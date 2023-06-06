from hub.dataload.nde import NDESourceUploader
from utils.pubtator import standardize_data
from utils.utils import check_schema


class DDEUploader(NDESourceUploader):
    name = "dde"
    __metadata__ = {
        "src_meta": {
            "url": "https://discovery.biothings.io/api/dataset/",
            "license_url": "https://creativecommons.org/licenses/by/4.0/",
            "license": "Creative Commons Attribution 4.0 International",
        }
    }

    @check_schema
    def load_data(self, data_folder):
        pubtator_docs = standardize_data(data_folder)
        for doc in pubtator_docs:
            yield doc
