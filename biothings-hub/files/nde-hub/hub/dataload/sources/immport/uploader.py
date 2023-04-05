from hub.dataload.nde import NDESourceUploader
from utils.pubtator import standardize_data
from utils.funding_helper import standardize_funding


class ImmPortUploader(NDESourceUploader):
    name = "immport"
    __metadata__ = {
        "src_meta": {
            "url": "https://www.immport.org/shared/home",
            "license_url": "https://docs.immport.org/home/agreement/"
        }
    }

    def load_data(self, data_folder):
        pubtator_docs = standardize_data(data_folder)
        funding_docs = standardize_funding(pubtator_docs)
        for doc in funding_docs:
            yield doc
