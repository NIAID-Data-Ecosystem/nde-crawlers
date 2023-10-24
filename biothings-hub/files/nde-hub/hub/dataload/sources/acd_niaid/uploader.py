from hub.dataload.nde import NDESourceUploader
from utils.funding_helper import standardize_funding
from utils.pmid_helper import load_pmid_ctfd
from utils.pubtator import standardize_data
from utils.utils import nde_upload_wrapper


class AccessClinicalDataUploader(NDESourceUploader):
    name = "acd_niaid"

    __metadata__ = {
        "merger": "merge_struct",
        "merger_kwargs": {"aslistofdict": "includedInDataCatalog", "include": ["includedInDataCatalog"]},
    }

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        docs = standardize_funding(docs)
        docs = standardize_data(docs)
        for doc in docs:
            yield doc
