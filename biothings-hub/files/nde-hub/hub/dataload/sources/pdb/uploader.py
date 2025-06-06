from hub.dataload.nde import NDESourceUploader
from utils.corrections import corrections
from utils.pubtator import standardize_data
from utils.utils import nde_upload_wrapper

# Example __metadata__ dictionary:
# <SOURCE_NAME> = https://api.data.niaid.nih.gov/v1/metadata
# __metadata__ = {
# "src_meta": {
# 'description': 'A short description of what the source offers, usually found on the source's about page',
# 'name': 'The full source name, Ex. Mendeley Data (not mendeley)',
# 'identifier':'includedInDataCatalog.name value',
# 'schema': 'A dict where the key is the source's metadata variable and the value is our transformation. Ex: {"summary":"description"},
# 'url': 'The source's URL homepage',
# }
# }


class PDB_Uploader(NDESourceUploader):
    name = "pdb"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = corrections(data_folder)
        for doc in docs:
            yield doc
