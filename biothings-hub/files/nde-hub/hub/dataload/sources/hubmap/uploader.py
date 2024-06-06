from hub.dataload.nde import NDESourceUploader
from utils.csv_helper import get_source_data
from utils.topic_category_helper import add_topic_category
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


class Hubmap_Uploader(NDESourceUploader):
    name = "hubmap"
    __metadata__ = {
        "src_meta": {
            "sourceInfo": {
                "description": "HuBMAP is part of a rich ecosystem of established and emerging atlasing programs supported by NIH and globally by other funding organizations, many of which are focused on specific organs or diseases. HuBMAP has connected with these programs to ensure data interoperability, avoid duplication of work, and leverage and synergize gained knowledge. The consortium has organized a number of events to bring together these communities to discuss topics of shared interest and is committed to improving coordination and collaboration among different programs. In addition, many of the HuBMAP PIs had been or are still actively participating in these efforts, helping with cross-pollination and advancing our global understanding. HuBMAP, as its name implies, was specifically initiated to resolve the challenge of building integrated, comprehensive, high-resolution spatial maps of human tissues and organs, which has resulted in HuBMAP providing leadership in the ecosystem around techniques for integrating disparate, multi-dimensional and multi-scale datasets, the development of a Common Coordinate Framework (CCF) for integrating data across many individuals, and the development and validation of these assays. To further increase interoperability, HuBMAP has adopted a number of standards and processes developed by other domain expert consortia, working and is actively involved in the knowledge exchange. The consortium sees itself as an integral part of the ecosystem, sharing its strengths and actively contributing to the community.",
                "identifier": "HuBMAP",
                "name": "HuBMAP",
                "schema": get_source_data(name),
                "url": "https://hubmapconsortium.org/",
            }
        }
    }

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = add_topic_category(data_folder, self.name)
        for doc in docs:
            yield doc
