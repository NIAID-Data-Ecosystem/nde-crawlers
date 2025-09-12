import biothings.hub.dataload.uploader as uploader
from hub.dataload.nde import NDESourceUploader


class New_NCBI_Geo_Uploader(uploader.BaseSourceUploader):
    name = "new_ncbi_geo"

    def load_data(self, data_folder):
        yield {"_id": "test_id", "name": "test_name"}

    @classmethod
    def get_mapping(cls):
        return NDESourceUploader.get_mapping()
