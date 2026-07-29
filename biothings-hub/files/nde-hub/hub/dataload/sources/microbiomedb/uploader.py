from hub.dataload.nde import NDESourceUploader
from utils.csv_helper import get_source_data


class Microbiomedb_Uploader(NDESourceUploader):
    name = "microbiomedb"

    __metadata__ = {
        "src_meta": {
            "sourceInfo": {
                "description": "MicrobiomeDB was developed as a discovery tool that empowers researchers to fully leverage their experimental metadata to construct queries that interrogate microbiome datasets.",
                "name": "MicrobiomeDB",
                "identifier": "MicrobiomeDB",
                "schema": get_source_data(name),
                "url": "https://beta.microbiomedb.org/mbio.beta/app/",
            }
        }
    }
