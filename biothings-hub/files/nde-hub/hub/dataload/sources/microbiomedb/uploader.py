from hub.dataload.nde import NDESourceUploader
from utils.csv_helper import get_source_data
from utils.pmid_helper import load_pmid_ctfd
from utils.utils import nde_upload_wrapper
from utils.pubtator import standardize_data

class Microbiomedb_Uploader(NDESourceUploader):
    name = "microbiomedb"

    __metadata__ = {
        "src_meta": {
            "sourceInfo": {
                "description": "MicrobiomeDB was developed as a discovery tool that empowers researchers to fully leverage their experimental metadata to construct queries that interrogate microbiome datasets.",
                "name": "MicrobiomeDB",
                "identifier": "MicrobiomeDB",
                "schema": get_source_data(name),
                "url": "https://beta.microbiomedb.org/mbio.beta/app/",
            }
        }
    }

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        pubtator_docs = standardize_data(docs)
        for doc in pubtator_docs:
            yield doc
