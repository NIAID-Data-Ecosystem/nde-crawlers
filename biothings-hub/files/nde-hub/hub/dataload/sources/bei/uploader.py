from hub.dataload.nde import NDESourceSampleUploader
from utils.extract import process_descriptions
from utils.pubtator import standardize_data
from utils.utils import nde_upload_wrapper


class BeiUploader(NDESourceSampleUploader):
    name = "bei"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = standardize_data(data_folder)
        docs = process_descriptions(docs)
        for doc in docs:
            yield doc
