from hub.dataload.nde import NDESourceUploader
from utils.pmid_helper import load_pmid_ctfd

class OmicsDIUploader(NDESourceUploader):
    name = "omicsdi"

    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        for doc in docs:
            yield doc