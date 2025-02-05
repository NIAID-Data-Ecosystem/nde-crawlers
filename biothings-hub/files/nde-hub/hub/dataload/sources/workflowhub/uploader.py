from hub.dataload.nde import NDESourceUploader
from utils.corrections import corrections
from utils.csv_helper import get_source_data
from utils.funding_helper import standardize_funding

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


class WorkflowHub_Uploader(NDESourceUploader):
    name = "workflowhub"
    __metadata__ = {
        "src_meta": {
            "sourceInfo": {
                "description": "WorkflowHub is a registry for describing, sharing and publishing scientific computational workflows. The registry supports any workflow in its native repository. WorkflowHub aims to facilitate discovery and re-use of workflows in an accessible and interoperable way. This is achieved through extensive use of open standards and tools, including Common Workflow Language (CWL), RO-Crate, BioSchemas and TRS, in accordance with the FAIR principles.",
                "identifier": "WorkflowHub",
                "name": "WorkflowHub",
                "schema": get_source_data(name),
                "url": "https://workflowhub.eu/",
            }
        }
    }
