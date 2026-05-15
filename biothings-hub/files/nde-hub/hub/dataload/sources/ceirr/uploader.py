from hub.dataload.nde import NDESourceSampleUploader
from utils.extract import process_descriptions
from utils.funding_helper import standardize_funding
from utils.pubtator import standardize_data
from utils.utils import nde_upload_wrapper


class CeirrUploader(NDESourceSampleUploader):
    name = "ceirr"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = standardize_funding(data_folder)
        docs = standardize_data(docs)
        docs = process_descriptions(docs)
        for doc in docs:
            yield doc
