from hub.dataload.nde import NDESourceUploader
from utils.corrections import load_documents
from utils.funding_helper import standardize_funding
from utils.pubtator import standardize_data
from utils.utils import nde_upload_wrapper


class DBAASP_Uploader(NDESourceUploader):
    name = "dbaasp"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_documents(data_folder)
        docs = standardize_data(docs)
        docs = standardize_funding(docs)
        for doc in docs:
            yield doc
