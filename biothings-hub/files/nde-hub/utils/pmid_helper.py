# Helper file to batch call pmids to get citations and funding
# Helpful links to documentation of Biopython Package for writing this file
# https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec162
# https://biopython.org/docs/1.76/api/Bio.Entrez.html
# https://www.nlm.nih.gov/bsd/mms/medlineelements.html
# https://dataguide.nlm.nih.gov/eutilities/utilities.html#efetch
import os
import orjson
import yaml
import time

from .date import add_date
from Bio import Entrez
from Bio import Medline
from datetime import datetime
from typing import Optional, List, Iterable, Dict
from itertools import islice
from config import GEO_API_KEY, GEO_EMAIL

import logging
logger = logging.getLogger('nde-logger')

def get_pub_date(date: str):
    """helper method to solve the problem transforming dates such as "2000 Spring" into isoformat dates

        https://www.nlm.nih.gov/bsd/mms/medlineelements.html#dp
        Returns:
            An isoformatted date depending on context:

            Seasons use metrological start
            Winter: December 1
            Spring: March 1
            Summer: June 1
            Fall: September 1

            Dates with Y/M/D-D take only the beginning day
            Dates with only Y/M or Y/M-M take the first day of that month
            Dates with only Y or Y-Y take first day of that year
            TODO: Not important since only one instance so far but fix edge case "2016 11-12"
    """

    months = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
    seasons = {"spring": " mar 1", "summer": " jun 1", "fall": " sep 1", "winter": " dec 1"}

    s_date = date.lower().split()
    date_len = len(s_date)
    # if length is 1 can either be year or year-year
    if date_len == 1:
        return datetime.strptime(s_date[0].split('-')[0], '%Y').isoformat()
    # if length is 2 can either be year season or year month or year month-month
    elif date_len == 2:
        if s_date[1][:3] in months:
            return datetime.strptime(s_date[0] + ' ' + s_date[1][:3], '%Y %b').isoformat()
        elif season := seasons.get(s_date[1]):
            return datetime.strptime(s_date[0] + season, '%Y %b %d').isoformat()
        else:
            logger.warning("Need to update isoformat transformation: %s", date)
            return None
    # if length is 3 should be year month day or year month day-day
    elif date_len == 3:
        return datetime.strptime(s_date[0] + ' ' + s_date[1] + ' ' + s_date[2].split('-')[0], '%Y %b %d').isoformat()
    # exception case there are quite a few entries with this case "2020 Jan - Feb"
    elif date_len == 4:
        if s_date[1] in months and s_date[3] in months and s_date[2] == "-":
            return datetime.strptime(s_date[0] + ' ' + s_date[1], '%Y %b').isoformat()
        else:
            logger.warning("Need to update isoformat transformation %s", date)
    else:
        logger.warning("Need to update isoformat transformation: %s", date)
        return None


# TODO add retry decorator to this function if getting too many imcompleteread errors
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

    # api query to parse citations
    handle = Entrez.efetch(db="pubmed", id=pmids, rettype="medline", retmode="text")

    records = Medline.parse(handle)
    # TODO: this can get an incompleteread error need implementation to rerun api query if this happens
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
    # TODO: this can get an incompleteread error need implementation to rerun api query if this happens
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

def load_pmid_ctfd(data_folder):
    """ Takes 1000 documents at a time and batch queries all of the pmids in the documents to improve runtime.
        Loads the citation and funding into the documents along with uploading the date field.
          Returns: A generator with the completed documents
    """

    # a temporary solution to make bigger batch api call instead of multiple smaller calls in crawler to improve runtime
    # TODO: figure out how to make a batch api call in crawler perferrably

    api_key = GEO_API_KEY
    email = GEO_EMAIL
    
    # if no access to config file comment out above and enter your own email
    # email = myemail@gmail.com

    with open(os.path.join(data_folder, 'data.ndjson'), 'rb') as f:
        while True:
            # pmid list for batch query
            pmid_list = []
            # docs to yield for each batch query
            doc_list = []

            # to make batch api query take the next 1000 docs and collect all the pmids
            next_n_lines = list(islice(f, 1000))
            if not next_n_lines:
                break
            for line in next_n_lines:
                doc = orjson.loads(line)
                doc_list.append(doc)
                if pmids := doc.get('pmids'):
                    pmid_list += [pmid.strip() for pmid in pmids.split(',')]

            # batch request
            eutils_info = batch_get_pmid_eutils(pmid_list, email, api_key)
            # throttle request rates, NCBI says up to 10 requests per second with API Key, 3/s without.
            if api_key:
                time.sleep(0.1)
            else:
                time.sleep(0.35)

            # add in the citation and funding to each doc in doc_list and yield
            for rec in doc_list:
                if pmids := rec.pop('pmids', None):
                    pmids = [pmid.strip() for pmid in pmids.split(',')]
                    # fixes issue where pmid numbers under 10 is read as 04 instead of 4
                    pmids = [ele.lstrip('0') for ele in pmids]
                    for pmid in pmids:
                        if not eutils_info.get(pmid):
                            logger.info('There is an issue with this pmid. PMID: %s, rec_id: %s', pmid, rec['_id'])
                        # this fixes the error where there is no pmid
                        # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41964
                        if eutils_info.get(pmid):
                            if citation := eutils_info[pmid].get('citation'):
                                if rec.get('citation'):
                                    rec['citation'].append(citation)
                                else:
                                    rec['citation'] = [citation]
                            if funding := eutils_info[pmid].get('funding'):
                                if rec.get('funding'):
                                    rec['funding'] += funding
                                else:
                                    rec['funding'] = funding
                # add the date tranformation before yielding
                rec = add_date(rec)
                yield rec

