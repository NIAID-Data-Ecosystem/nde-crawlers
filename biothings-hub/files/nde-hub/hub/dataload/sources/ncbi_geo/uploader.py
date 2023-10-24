from hub.dataload.nde import NDESourceUploader
from utils.funding_helper import standardize_funding
from utils.pmid_helper import load_pmid_ctfd
from utils.pubtator import standardize_data
from utils.utils import nde_upload_wrapper


class NCBI_Geo_Uploader(NDESourceUploader):
    name = "ncbi_geo"

    __metadata__ = {
        "src_meta": {
            "url": "https://www.ncbi.nlm.nih.gov/geo/browse/",
            "license_url": "https://www.ncbi.nlm.nih.gov/home/about/policies/",
        },
        "merger": "merge_struct",
        "merger_keywords": {"aslistofdict": "includedInDataCatalog", "include": ["includedInDataCatalog"]},
    }

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        docs = standardize_funding(docs)
        pubtator_docs = standardize_data(docs)
        for doc in pubtator_docs:
            yield doc
