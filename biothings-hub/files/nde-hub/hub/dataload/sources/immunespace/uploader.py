from hub.dataload.nde import NDESourceUploader
from utils import iter_ndjson, nde_upload_wrapper


class ImmunespaceUploader(NDESourceUploader):
    # TODO metadata description
    __metadata__ = {
        "merger": "merge_struct",
        "merger_kwargs": {"aslistofdict": "includedInDataCatalog", "include": ["includedInDataCatalog"]},
    }

    name = "immunespace"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        yield from (doc for doc in iter_ndjson(data_folder) if doc.get("@type") == "Dataset")
