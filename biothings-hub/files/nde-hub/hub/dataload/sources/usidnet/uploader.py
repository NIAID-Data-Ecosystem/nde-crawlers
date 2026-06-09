from hub.dataload.nde import NDESourceSampleUploader
from utils.pubtator import standardize_data
from utils.utils import nde_upload_wrapper


class UsidnetUploader(NDESourceSampleUploader):
    name = "usidnet"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = standardize_data(data_folder)
        for doc in docs:
            yield doc
