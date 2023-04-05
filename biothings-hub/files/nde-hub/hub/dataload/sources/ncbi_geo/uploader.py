from hub.dataload.nde import NDESourceUploader
from utils.pmid_helper import load_pmid_ctfd
from utils.utils import check_schema


class NCBI_Geo_Uploader(NDESourceUploader):
    name = "ncbi_geo"

    __metadata__ = {
        "src_meta": {
            "url": "https://www.ncbi.nlm.nih.gov/geo/browse/",
            "license_url": "https://www.ncbi.nlm.nih.gov/home/about/policies/"
        }
    }

    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        for doc in docs:
            # check schema
            check_schema(doc)
            yield doc
