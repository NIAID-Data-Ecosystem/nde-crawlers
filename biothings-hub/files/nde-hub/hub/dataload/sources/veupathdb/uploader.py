from hub.dataload.nde import NDESourceUploader
from utils.pmid_helper import load_pmid_ctfd
from utils.utils import check_schema
from utils.pubtator import standardize_data


class VEuPathDB_Uploader(NDESourceUploader):
    name = "veupathdb"

    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)

        pubtator_docs = standardize_data(docs)
        for doc in pubtator_docs:
            # check schema
            check_schema(doc)
            yield doc
