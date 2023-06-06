from hub.dataload.nde import NDESourceUploader
from utils.pmid_helper import load_pmid_ctfd
from utils.utils import check_schema


class OmicsDIUploader(NDESourceUploader):
    name = "omicsdi"
    __metadata__ = {
        "src_meta": {"url": "https://www.omicsdi.org/search", "license_url": "https://www.ebi.ac.uk/licencing"}
    }

    @check_schema
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        for doc in docs:
            yield doc
