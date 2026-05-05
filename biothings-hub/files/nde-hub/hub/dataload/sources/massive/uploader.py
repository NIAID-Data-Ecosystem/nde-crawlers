from hub.dataload.nde import NDESourceUploader
from utils.extract import process_descriptions
from utils.measurement_technique_helper import process_measurement_technique
from utils.pubtator import standardize_data
from utils.utils import nde_upload_wrapper


class MassiveUploader(NDESourceUploader):
    name = "massive"

    __metadata__ = {
        "merger": "merge_struct",
        "merger_kwargs": {"aslistofdict": "includedInDataCatalog", "include": ["includedInDataCatalog"]},
    }

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = standardize_data(data_folder)
        docs = process_descriptions(docs)
        docs = process_measurement_technique(docs, self.name)
        for doc in docs:
            yield doc
