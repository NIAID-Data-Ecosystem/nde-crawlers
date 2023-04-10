import os
import orjson
from utils.utils import check_schema
from hub.dataload.nde import NDESourceUploader


class ZenodoUploader(NDESourceUploader):
    name = "zenodo"

    @check_schema
    def load_data(self, data_folder):
        with open(os.path.join(data_folder, 'data.ndjson'), 'rb') as f:
            for line in f:
                doc = orjson.loads(line)
                yield doc