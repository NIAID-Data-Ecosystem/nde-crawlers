from hub.dataload.nde import NDESourceUploader
from utils.corrections import corrections
from utils.csv_helper import get_source_data
from utils.extract import process_descriptions
from utils.funding_helper import standardize_funding
from utils.topic_category_helper import add_topic_category
from utils.utils import nde_upload_wrapper


class DataverseUploader(NDESourceUploader):
    name = "dataverse"
    __metadata__ = {
        "src_meta": {
            "sourceInfo": {
                "description": "The Harvard Dataverse Repository is a free data repository open to all researchers from any discipline, both inside and outside of the Harvard community, where you can share, archive, cite, access, and explore research data. Each individual Dataverse collection is a customizable collection of datasets (or a virtual repository) for organizing, managing, and showcasing datasets.",
                "identifier": "Harvard Dataverse",
                "name": "Harvard Dataverse",
                "schema": get_source_data(name),
                "url": "https://dataverse.harvard.edu/",
            }
        }
    }

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = standardize_funding(data_folder)
        docs = process_descriptions(docs)
        docs = corrections(docs)
        docs = add_topic_category(docs, self.name)
        for doc in docs:
            yield doc
