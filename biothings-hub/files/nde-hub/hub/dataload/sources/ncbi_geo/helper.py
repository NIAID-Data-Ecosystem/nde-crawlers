import time

from Bio import Entrez
from Bio import Medline
from datetime import datetime
from typing import Optional, List, Iterable, Dict

import logging

logger = logging.getLogger('nde-logger')

"""helper method to solve the problem transforming dates such as "2000 Spring" into isoformat dates

https://www.nlm.nih.gov/bsd/mms/medlineelements.html#dp
Returns:
    An isoformatted date depending on context:

    Seasons use metrological start
    Winter: December 1
    Spring: March 1
    Summer: June 1
    Fall: September 1

    Dates with only Y/M or Y/M-M take the first day of that month
    Dates with only Y take first day of that year
"""


def get_pub_date(date: str):
    months = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
    seasons = {"spring": " mar 1", "summer": " jun 1", "fall": " sep 1", "winter": " dec 1"}

    s_date = date.lower().split()
    date_len = len(s_date)
    # if length is 1 just the year
    if date_len == 1:
        return datetime.strptime(date, '%Y').isoformat()
    # if length is 2 can either be year season or year month
    elif date_len == 2:
        if s_date[1][:3] in months:
            return datetime.strptime(s_date[0] + ' ' + s_date[1][:3], '%Y %b').isoformat()
        elif season := seasons.get(s_date[1]):
            return datetime.strptime(s_date[0] + season, '%Y %b %d').isoformat()
        else:
            logger.warning("Need to update isoformat transformation: %s", date)
            return None
    # if length is 3 should be year month day
    elif date_len == 3:
        return datetime.strptime(date, '%Y %b %d').isoformat()
    else:
        logger.warning("Need to update isoformat transformation: %s", date)
        return None


def batch_get_pmid_eutils(pmids: Iterable[str], email: str, api_key: Optional[str] = None) -> Dict:
    """Use pmid to retrieve both citation and funding info in batch
    :param pmids: A list of PubMed PMIDs
    :param api_key: API Key from NCBI to access E-utilities
    :return: A dictionary containing the pmids which hold citations and funding.
    """
    # probably dont need this line. Using python package, should work both ways.
    # if pmids is str:
    #     warnings.warn(f"Got str:{pmids} as parameter, expecting an Iterable of str", RuntimeWarning)

    # set up Entrez variables. Email is required.
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    ct_fd = {}

    handle = Entrez.efetch(db="pubmed", id=pmids, rettype="medline", retmode="text")

    records = Medline.parse(handle)
    records = list(records)
    for record in records:
        citation = {}
        # rename keys
        if name := record.get('TI'):
            citation['name'] = name
        if pmid := record.get('PMID'):
            citation['pmid'] = pmid
            citation['identifier'] = 'PMID:' + pmid
            citation['url'] = 'https://pubmed.ncbi.nlm.nih.gov/' + pmid + '/'
        if journal_name := record.get('JT'):
            citation['journalName'] = journal_name
        if date_published := record.get('DP'):
            if date := get_pub_date(date_published):
                citation['datePublished'] = date

        # make an empty list if there is some kind of author
        if record.get('AU') or record.get('CN'):
            citation['author'] = []
        if authors := record.get('AU'):
            for author in authors:
                citation['author'].append({'@type': 'Person', 'name': author})
        if corp_authors := record.get('CN'):
            for corp_author in corp_authors:
                citation['author'].append({'@type': 'Organization', 'name': corp_author})
        # put citation in dictionary
        ct_fd[pmid] = {'citation': citation}

    # throttle request rates, NCBI says up to 10 requests per second with API Key, 3/s without.
    if api_key:
        time.sleep(0.1)
    else:
        time.sleep(0.35)

    # get the funding using xml file because of problems parsing the medline file
    # https://www.nlm.nih.gov/bsd/mms/medlineelements.html#gr
    handle = Entrez.efetch(db="pubmed", id=pmids, retmode="xml")

    # Have to use Entrez.read() instead of Entrez.parse(). The problem is discussed here: https://github.com/biopython/biopython/issues/1027
    records = Entrez.read(handle)
    records = records['PubmedArticle']

    funding = []
    for record in records:
        if grants := record['MedlineCitation']['Article'].get('GrantList'):
            for grant in grants:
                fund = {}
                if grant_id := grant.get('GrantID'):
                    fund['identifier'] = grant_id
                if agency := grant.get('Agency'):
                    fund['funder'] = {'@type': 'Organization', 'name': agency}
                funding.append(fund)
        if pmid := record['MedlineCitation'].get('PMID'):
            if funding:
                ct_fd[pmid]['funding'] = funding
            funding = []

    return ct_fd
