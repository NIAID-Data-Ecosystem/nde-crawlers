from hub.dataload.nde import NDESourceUploader
from utils.funding_helper import standardize_funding
from utils.pubtator import standardize_data
from utils.topic_category_helper import add_topic_category
from utils.utils import nde_upload_wrapper


class ImmPortUploader(NDESourceUploader):
    name = "immport"
    __metadata__ = {
        "src_meta": {
            "url": "https://www.immport.org/shared/home",
            "license_url": "https://docs.immport.org/home/agreement/",
        }
    }

    @nde_upload_wrapper
    def load_data(self, data_folder):
        pubtator_docs = standardize_data(data_folder)
        funding_docs = standardize_funding(pubtator_docs)
        docs = add_topic_category(funding_docs, self.name)
        for doc in docs:
            yield doc
