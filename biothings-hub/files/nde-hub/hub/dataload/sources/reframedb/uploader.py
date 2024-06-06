from hub.dataload.nde import NDESourceUploader
from utils.funding_helper import standardize_funding
from utils.pmid_helper import load_pmid_ctfd
from utils.pubtator import standardize_data
from utils.topic_category_helper import add_topic_category
from utils.utils import nde_upload_wrapper


class ReframedbUploader(NDESourceUploader):
    name = "reframedb"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        docs = standardize_funding(docs)
        pubtator_docs = standardize_data(docs)
        topic_category_docs = add_topic_category(pubtator_docs, self.name)
        for doc in topic_category_docs:
            yield doc
