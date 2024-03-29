from hub.dataload.nde import NDESourceUploader
from utils.pmid_helper import load_pmid_ctfd
from utils.pubtator import standardize_data
from utils.utils import nde_upload_wrapper
from utils.funding_helper import standardize_funding


class ReframedbUploader(NDESourceUploader):
    name = "reframedb"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        docs = standardize_funding(docs)
        pubtator_docs = standardize_data(docs)
        for doc in pubtator_docs:
            yield doc
