import os
import orjson
import yaml
import time

from .helper import batch_get_pmid_eutils
from hub.dataload.nde import NDESourceUploader
from itertools import islice
from config import GEO_API_KEY, GEO_EMAIL

class NCBI_Geo_Uploader(NDESourceUploader):
    name = "ncbi_geo"

    def load_data(self, data_folder):
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
                        for pmid in pmids:
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
                    yield rec

