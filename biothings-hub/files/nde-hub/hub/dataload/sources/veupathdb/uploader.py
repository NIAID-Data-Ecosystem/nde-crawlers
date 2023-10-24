from hub.dataload.nde import NDESourceUploader
from utils.funding_helper import standardize_funding
from utils.pmid_helper import load_pmid_ctfd
from utils.pubtator import standardize_data
from utils.utils import nde_upload_wrapper


class VEuPathDB_Uploader(NDESourceUploader):
    name = "veupathdb"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        docs = standardize_funding(docs)
        docs = standardize_data(docs)
        for doc in docs:
            yield doc
