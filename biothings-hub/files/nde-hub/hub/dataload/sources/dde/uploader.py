from hub.dataload.nde import NDESourceUploader
from utils.in_defined_term_set import handle_dde_docs
from utils.pubtator import standardize_data
from utils.utils import nde_upload_wrapper


class DDEUploader(NDESourceUploader):
    name = "dde"
    __metadata__ = {
        "src_meta": {
            "url": "https://discovery.biothings.io/api/dataset/",
            "license_url": "https://creativecommons.org/licenses/by/4.0/",
            "license": "Creative Commons Attribution 4.0 International",
        }
    }

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = handle_dde_docs(data_folder)
        pubtator_docs = standardize_data(docs)
        for doc in pubtator_docs:
            yield doc
