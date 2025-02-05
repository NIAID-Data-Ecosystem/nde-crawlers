from hub.dataload.nde import NDESourceUploader
from utils.corrections import corrections
from utils.csv_helper import get_source_data
from utils.funding_helper import standardize_funding
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


class Biocontainers_Uploader(NDESourceUploader):
    name = "biocontainers"
    __metadata__ = {
        "src_meta": {
            "sourceInfo": {
                "description": "BioContainers is a community-driven project that provides the infrastructure and basic guidelines to create, manage and distribute bioinformatics containers with a special focus on omics fields such as proteomics, genomics, transcriptomics and metabolomics. The currently available BioContainers containers facilitate the usage, and reproducibility of software and algorithms. They can be integrated into more comprehensive bioinformatics pipelines and different architectures (local desktop, Cloud environments or HPC clusters). BioContainers is based on the popular frameworks Conda, Docker and Singularity.",
                "identifier": "BioContainers",
                "name": "BioContainers",
                "schema": get_source_data(name),
                "url": "https://biocontainers.pro/",
            }
        }
    }

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        docs = standardize_funding(docs)
        docs = corrections(docs)
        for doc in docs:
            yield doc
