from hub.dataload.nde import NDESourceUploader
from utils.csv_helper import get_source_data
from utils.pmid_helper import load_pmid_ctfd
from utils.utils import check_schema

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


class Qiita_Uploader(NDESourceUploader):
    name = "qiita"
    __metadata__ = {
        "src_meta": {
            "sourceInfo": {
                "description": "Qiita(canonically pronounced cheetah) is an entirely open-source microbial study management platform. It allows users to keep track of multiple studies with multiple 'omics data. Additionally, Qiita is capable of supporting multiple analytical pipelines through a 3rd-party plugin system, allowing the user to have a single entry point for all of their analyses. Qiita provides database and compute resources to the global community, alleviating the technical burdens that are typically limiting for researchers studying microbial ecology(e.g. familiarity with the command line or access to compute power).Qiita's platform allows for quick reanalysis of the datasets that have been deposited using the latest analytical technologies. This means that Qiita's internal datasets are living data that is periodically re-annotated according to current best practices.",
                "identifier": "Qiita",
                "name": "Qiita",
                "schema": get_source_data(name),
                "url": "https://qiita.ucsd.edu/"
            }
        }
    }

    @check_schema
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        for doc in docs:
            yield doc
