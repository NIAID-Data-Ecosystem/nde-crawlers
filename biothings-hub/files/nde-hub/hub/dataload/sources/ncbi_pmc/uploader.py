from hub.dataload.nde import NDESourceUploader
from utils.pmid_helper import load_pmid_ctfd
from utils.utils import nde_upload_wrapper
from utils.pubtator import standardize_data


class NCBI_PMC_Uploader(NDESourceUploader):
    name = "ncbi_pmc"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        pubtator_docs = standardize_data(docs)
        for doc in pubtator_docs:
            yield doc
