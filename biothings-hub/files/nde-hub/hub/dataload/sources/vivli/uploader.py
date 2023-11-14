from hub.dataload.nde import NDESourceUploader
from utils.clinical_trails_helper import load_ct_wrapper
from utils.funding_helper import standardize_funding
from utils.pubtator import standardize_data
from utils.utils import nde_upload_wrapper


class VivliUploader(NDESourceUploader):
    name = "vivli"
    __metadata__ = {
        "src_meta": {
            "url": "https://search.vivli.org/",
            "license_url": "https://vivli.org/resources/vivli-data-use-agreement/",
        }
    }

    @nde_upload_wrapper
    @load_ct_wrapper
    def load_data(self, data_folder):
        docs = standardize_data(data_folder)
        docs = standardize_funding(docs)
        for doc in docs:
            yield doc
