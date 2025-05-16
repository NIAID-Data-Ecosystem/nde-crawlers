from hub.dataload.nde import NDESourceUploader
from utils.corrections import corrections
from utils.extract import process_descriptions
from utils.pubtator import standardize_data
from utils.utils import nde_upload_wrapper


class ImmunespaceUploader(NDESourceUploader):
    # TODO metadata description
    __metadata__ = {
        "merger": "merge_struct",
        "merger_kwargs": {"aslistofdict": "includedInDataCatalog", "include": ["includedInDataCatalog"]},
    }

    name = "immunespace"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = standardize_data(data_folder)
        docs = corrections(docs)

        for doc in docs:
            yield doc
