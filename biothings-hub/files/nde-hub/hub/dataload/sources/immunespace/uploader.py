from hub.dataload.nde import NDESourceUploader
from utils.corrections import corrections
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
        docs = corrections(data_folder)
        for doc in docs:
            yield doc
