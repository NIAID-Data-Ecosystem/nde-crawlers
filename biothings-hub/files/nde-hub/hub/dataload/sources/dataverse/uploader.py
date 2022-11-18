from hub.dataload.nde import NDESourceUploader
from utils.pmid_helper import load_pmid_ctfd
from utils.csv_helper import get_source_data

class DataverseUploader(NDESourceUploader):
    name = "dataverse"
    __metadata__ = {
        "src_meta": {
            "sourceInfo": {
                "description": "The Harvard Dataverse Repository is a free data repository open to all researchers from any discipline, both inside and outside of the Harvard community, where you can share, archive, cite, access, and explore research data. Each individual Dataverse collection is a customizable collection of datasets (or a virtual repository) for organizing, managing, and showcasing datasets.",
                "identifier": "Dataverse",
                "name": "Dataverse",
                "schema": get_source_data(name),
                "url": "https://dataverse.harvard.edu/"
            }
        }
    }

    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        for doc in docs:
            yield doc
