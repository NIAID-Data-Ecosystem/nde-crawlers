from hub.dataload.nde import NDESourceUploader
from utils.corrections import load_documents
from utils.funding_helper import standardize_funding
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


class EMDB_Uploader(NDESourceUploader):
    name = "emdb"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_documents(data_folder)
        docs = standardize_data(docs)
        docs = standardize_funding(docs)
        for doc in docs:
            yield doc
