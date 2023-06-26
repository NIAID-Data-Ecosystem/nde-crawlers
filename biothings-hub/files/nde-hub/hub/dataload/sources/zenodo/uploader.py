import os

import orjson
from hub.dataload.nde import NDESourceUploader
from utils.utils import nde_upload_wrapper


class ZenodoUploader(NDESourceUploader):
    name = "zenodo"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        with open(os.path.join(data_folder, "data.ndjson"), "rb") as f:
            for line in f:
                doc = orjson.loads(line)
                yield doc
