from hub.dataload.nde import NDESourceUploader
from utils.pmid_helper import load_pmid_ctfd
from utils.utils import nde_upload_wrapper
from utils.funding_helper import standardize_funding


class NCBI_PMC_Uploader(NDESourceUploader):
    name = "ncbi_pmc"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        docs = standardize_funding(docs)
        for doc in docs:
            yield doc
