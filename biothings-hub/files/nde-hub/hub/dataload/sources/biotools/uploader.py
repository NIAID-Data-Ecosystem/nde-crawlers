from hub.dataload.nde import NDESourceUploader
from utils.extract import process_descriptions
from utils.funding_helper import standardize_funding
from utils.utils import nde_upload_wrapper


class Biotools_Uploader(NDESourceUploader):
    name = "biotools"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = standardize_funding(data_folder)
        docs = process_descriptions(docs)
        for doc in docs:
            yield doc
