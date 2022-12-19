import json
import datetime
import requests

import logging

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(name)s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger('nde-logger')


def parse():


    url = 'https://beta.microbiomedb.org/mbio.beta/service/record-types/dataset/searches/AllDatasets/reports/standard?' \
              'reportConfig={"attributes":[acknowledgement,pmids,study_categories,description,type,institution,' \
              'card_headline,contact,project_availability,study_access,email,is_species_scope,summary,dataset_sha1_digest,' \
              'card_points,dataset_id,build_number_introduced,short_display_name,display_name,dataset_name,primary_key,version,' \
              'pmids_download,bulk_download_url,eda_study_id,short_attribution,dataset_display_name,is_public,' \
              'eupath_release,caveat],"tables":[Contacts,DownloadVersion,Publications,AssociatedDatasets,HyperLinks],' \
              '"attributeFormat":"text"}'

    logger.info("Making request...")
    request = requests.get(url)
    logger.info("Request made. HTTP STATUS: %d", request.status_code)
    records = request.json()
    logger.info("Parsing records...")
    for count, record in enumerate(records['records'], start=1):
        if count % 1000 == 0:
            logger.info("Parsed %d records", count)

        output = {
            'includedInDataCatalog': {
                '@type': 'Dataset',
                'name': 'MicrobiomeDB',
                'url': 'https://beta.microbiomedb.org/mbio.beta/app/',
                'versionDate': datetime.date.today().isoformat()
            },
            '@type': 'Dataset',
            '_id': "MICROBIOME_" + record['attributes'].get('dataset_id')
        }

        # set up keywords list
        k = []
        if keywords := record['attributes'].get('study_categories'):
            keywords = json.loads(keywords)
            k = keywords + k

        # type should only be a string
        if keyword := record['attributes'].get('type'):
            k.append(keyword)

        if k:
            output['keywords'] = k

        if description := record['attributes'].get('description'):
            output['description'] = description

        # set up contributor field
        contributor = {}
        if contributor_aff := record['attributes'].get('institution'):
            contributor['affiliation'] = contributor_aff

        if contributor_name := record['attributes'].get('name'):
            contributor['name'] = contributor_name

        if contributor_email := record['attributes'].get('email'):
            contributor['email'] = contributor_email

        if contributor:
            output['contributor'] = contributor

        if coa := record['attributes'].get('study_access'):
            output['conditionsOfAccess'] = coa

        if abstract := record['attributes'].get('summary'):
            output['abstract'] = abstract

        if description := record['attributes'].get('card_points'):
            description = json.loads(description)
            output['description'] = description

        if alternate_name := record['attributes'].get('dataset_name'):
            output['alternateName'] = alternate_name

        if pmids := record['attributes'].get('pmids_download'):
            output['pmids'] = pmids

        # if distribution_content_url := record['attributes'].get('bulk_url_download'):
        #     output['distribution'] = {'contentUrl': distribution_content_url}

        # set up author field
        if authors := record['tables'].get('Contacts'):
            a_list = []
            for author in authors:
                a = {}
                if name := author.get('contact_name'):
                    a['name'] = name
                if affiliation := author.get('affiliation'):
                    a['affiliation'] = {'name': affiliation}
                if a:
                    a_list.append(a)
            if a_list:
                output['author'] = a_list

        # set up distribution field
        if distributions := record['tables'].get('DownloadVersion'):
            d_list = []
            for distribution in distributions:
                d = {}
                if date := distribution.get('release_date'):
                    date = datetime.datetime.strptime(date, '%Y-%b-%d').date().isoformat()
                    d['datePublished'] = date
                if identifier := distribution.get('dataset_id'):
                    d['identifier'] = identifier
                if name := distribution.get('dataset_name'):
                    d['name'] = name
                if d:
                    d_list.append(d)
            if d_list:
                output['distribution'] = d_list

        # get both isPartOf and isRelatedTo
        if ads := record['tables'].get('AssociatedDatasets'):
            # isPartOf
            ipo_list = []
            # isRelatedTo
            irt_list = []
            for ad in ads:
                ipo = {}
                if name := ad.get('project'):
                    ipo['name'] = name
                if alternate_name := ad.get('webapp_name'):
                    ipo['alternateName'] = alternate_name
                if ipo:
                    ipo_list.append(ipo)

                irt = {}
                if cite_url := ad.get('pubmed link'):
                    irt['citation'] = {'url': cite_url}

                if (cite_pmid := ad.get('pmid')) and cite_url:
                    irt['citation']['pmid'] = cite_pmid
                elif cite_pmid := ad.get('pmid'):
                    irt['citation'] = {'pmid': cite_pmid}

                if url := ad.get('associated dataset'):
                    irt['url'] = url
                if name := ad.get('display_name'):
                    irt['name'] = name
                if identifier := ad.get('associated_dataset_id'):
                    irt['identifier'] = identifier
                if irt:
                    irt_list.append(irt)
            if ipo_list:
                output['isPartOf'] = ipo_list
            if irt_list:
                output['isRelatedTo'] = irt_list

        # grabs urls
        if links := record['tables'].get('HyperLinks'):
            link_list = []
            for link in links:
                if url := link.get('url'):
                    link_list.append(url)
            if link_list:
                output['url'] = link_list

        yield output

    logger.info("Finished Parsing. Total records: %d", count)
