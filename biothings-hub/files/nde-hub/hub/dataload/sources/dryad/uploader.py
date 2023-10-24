from hub.dataload.nde import NDESourceUploader
from utils.funding_helper import standardize_funding
from utils.utils import nde_upload_wrapper


class DryadUploader(NDESourceUploader):
    name = "dryad"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = standardize_funding(data_folder)
        for doc in docs:
            yield doc
