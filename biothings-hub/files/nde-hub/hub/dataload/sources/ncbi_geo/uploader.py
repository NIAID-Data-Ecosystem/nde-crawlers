from hub.dataload.nde import NDESourceUploader
from utils.pmid_helper import load_pmid_ctfd

class NCBI_Geo_Uploader(NDESourceUploader):
    name = "ncbi_geo"

    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        for doc in docs:
            yield doc
            