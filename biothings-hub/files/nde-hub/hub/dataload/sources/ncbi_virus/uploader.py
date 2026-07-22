from hub.dataload.nde import NDESourceUploader
from utils.corrections import load_documents
from utils.pubtator import standardize_data
from utils.utils import nde_upload_wrapper


class NCBI_VIRUS_Uploader(NDESourceUploader):
    name = "ncbi_virus"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_documents(data_folder)
        docs = standardize_data(docs)
        for doc in docs:
            yield doc
