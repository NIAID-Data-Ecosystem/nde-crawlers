from hub.dataload.nde import NDESourceUploader
from utils.corrections import corrections
from utils.funding_helper import standardize_funding
from utils.pmid_helper import load_pmid_ctfd
from utils.pubtator import standardize_data
from utils.topic_category_helper import add_topic_category
from utils.utils import nde_upload_wrapper


class dbGaP_Uploader(NDESourceUploader):
    name = "dbgap"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        docs = standardize_funding(docs)
        docs = standardize_data(docs)
        docs = add_topic_category(docs)
        docs = corrections(docs)
        for doc in docs:
            yield doc
