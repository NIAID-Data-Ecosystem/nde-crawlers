from hub.dataload.nde import NDESourceUploader
from utils.csv_helper import get_source_data
from utils.funding_helper import standardize_funding
from utils.utils import nde_upload_wrapper

class Biotools_Uploader(NDESourceUploader):
    name = "biotools"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = standardize_funding(data_folder)
        for doc in docs:
            yield doc
