from hub.dataload.nde import NDESourceUploader
from utils.csv_helper import get_source_data
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


class HCA_Uploader(NDESourceUploader):
    name = "hca"
    __metadata__ = {
        "src_meta": {
            "sourceInfo": {
                "description": "Thanks to new single cell genomics and spatial imaging technologies developed since the late 2000s and early 2010s, it is now possible to measure gene expression profiles in individual cells. These large scale data can be used with machine learning algorithms to decipher how the cells differ from and interact with their neighbors, and how they form and function in the tissue. This now allows scientists to identify and understand cell types in unprecedented detail, resolution and breadth. The Human Cell Atlas (HCA) is an international group of researchers using a combination of these new technologies to create cellular reference maps with the position, function and characteristics of every cell type in the human body.",
                "identifier": "Human Cell Atlas",
                "name": "Human Cell Atlas",
                "schema": get_source_data(name),
                "url": "https://www.humancellatlas.org/",
            }
        }
    }

    @nde_upload_wrapper
    def load_data(self, data_folder):
        pubtator_docs = standardize_data(data_folder)
        for doc in pubtator_docs:
            yield doc
