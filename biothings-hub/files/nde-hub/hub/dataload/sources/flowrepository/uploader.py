from hub.dataload.nde import NDESourceUploader
from utils.csv_helper import get_source_data
from utils.pmid_helper import load_pmid_ctfd
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


class FlowRepositoryUploader(NDESourceUploader):
    name = "flowrepository"
    __metadata__ = {
        "src_meta": {
            "sourceInfo": {
                "description": "FlowRepository is a database of flow cytometry experiments where you can query and download data collected and annotated according to the MIFlowCyt standard. It is primarily used as a data deposition place for experimental findings published in peer-reviewed journals in the flow cytometry field.",
                "identifier": "FlowRepository",
                "name": "FlowRepository",
                "schema": get_source_data(name),
                "url": "http://flowrepository.org/",
            }
        }
    }

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        for doc in docs:
            yield doc
