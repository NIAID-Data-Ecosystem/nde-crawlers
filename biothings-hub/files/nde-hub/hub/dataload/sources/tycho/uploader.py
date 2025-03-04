from hub.dataload.nde import NDESourceUploader
from utils.funding_helper import standardize_funding
from utils.pubtator import standardize_data
from utils.utils import nde_upload_wrapper


class TychoUploader(NDESourceUploader):
    name = "tycho"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = standardize_funding(data_folder)
        docs = standardize_data(docs)
        for doc in docs:
            yield doc
