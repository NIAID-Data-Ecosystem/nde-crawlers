from hub.dataload.nde import NDESourceUploader
from utils.corrections import corrections
from utils.pubtator import standardize_data
from utils.pmid_helper import load_pmid_ctfd
from utils.topic_category_helper import add_topic_category
from utils.utils import nde_upload_wrapper

class NDExUploader(NDESourceUploader):
    name = "ndex"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        docs = standardize_data(docs)
        docs = corrections(docs)
        docs = add_topic_category(docs, self.name)
        for doc in docs:
            yield doc
