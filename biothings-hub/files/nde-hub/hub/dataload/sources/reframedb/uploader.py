from hub.dataload.nde import NDESourceUploader
from utils.pmid_helper import load_pmid_ctfd
from utils.utils import check_schema

class ReframedbUploader(NDESourceUploader):
    name = "reframedb"

    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        for doc in docs:
            # check schema
            check_schema(doc)
            yield doc
