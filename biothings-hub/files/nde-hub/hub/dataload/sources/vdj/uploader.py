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


class VDJ_Uploader(NDESourceUploader):
    name = "vdj"
    __metadata__ = {
        "src_meta": {
            "sourceInfo": {
                "description": "VDJServer is a free, scalable resource for performing immune repertoire analysis and sharing data. VDJServer Community Data Portal is part of the AIRR Data Commons. Funded by a National Institute of Allergy and Infectious Diseases research grant (#1R01A1097403), the VDJServer project is led by The University of Texas Southwestern (UTSW) Medical Center in collaboration with the J. Craig Venter Institute and Yale University. The Texas Advanced Computing Center (TACC) at The University of Texas at Austin leads the cyberinfrastructure implementation, including the high performance computing (HPC) systems, storage, and software solutions.",
                "identifier": "VDJServer",
                "name": "vdj",
                "schema": get_source_data(name),
                "url": "https://vdj-staging.tacc.utexas.edu/community/"
            }
        }
    }

    @check_schema
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        for doc in docs:
            yield doc
