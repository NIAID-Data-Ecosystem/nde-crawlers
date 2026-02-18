from hub.dataload.nde import NDESourceUploader
from utils.extract import process_descriptions
from utils.funding_helper import standardize_funding
from utils.topic_category_helper import add_topic_category
from utils.utils import nde_upload_wrapper


class DryadUploader(NDESourceUploader):
    name = "dryad"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = standardize_funding(data_folder)
        docs = process_descriptions(docs)
        docs = add_topic_category(docs, self.name)
        for doc in docs:
            yield doc
