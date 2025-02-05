from hub.dataload.nde import NDESourceUploader
from utils.corrections import corrections
from utils.extract import process_descriptions
from utils.utils import nde_upload_wrapper


class MalariaGenUploader(NDESourceUploader):
    name = "malariagen"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = corrections(data_folder)
        docs = process_descriptions(docs)
        for doc in docs:
            yield doc
