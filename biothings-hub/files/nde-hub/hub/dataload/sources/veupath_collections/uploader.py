import os

import orjson
from hub.dataload.nde import NDESourceUploader


class VEuPathCollections_Uploader(NDESourceUploader):
    # TODO metadata description
    __metadata__ = {
        "merger": "merge_struct",
        "merger_kwargs": {"aslistofdict": "includedInDataCatalog", "include": ["includedInDataCatalog"]},
    }

    name = "veupath_collections"

    def load_data(self, data_folder):
        with open(os.path.join(data_folder, "data.ndjson"), "rb") as f:
            for line in f:
                doc = orjson.loads(line)
                yield doc
