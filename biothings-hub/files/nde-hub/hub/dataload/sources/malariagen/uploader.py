from hub.dataload.nde import NDESourceUploader
from utils.extract import process_descriptions
from utils.utils import nde_upload_wrapper


class MalariaGenUploader(NDESourceUploader):
    name = "malariagen"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = process_descriptions(data_folder)
        for doc in docs:
            yield doc
