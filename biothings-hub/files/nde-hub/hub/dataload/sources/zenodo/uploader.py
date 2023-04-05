import os
import orjson
from utils.date import add_date
from hub.dataload.nde import NDESourceUploader


class ZenodoUploader(NDESourceUploader):
    name = "zenodo"

    def load_data(self, data_folder):
        with open(os.path.join(data_folder, 'data.ndjson'), 'rb') as f:
            for line in f:
                doc = orjson.loads(line)
                # add date transformation here and have the most recent date of all of the dates
                doc = add_date(doc)
                yield doc