from hub.dataload.nde import NDESourceUploader
from utils.csv_helper import get_source_data
from utils.funding_helper import standardize_funding
from utils.pmid_helper import load_pmid_ctfd
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


class dashUploader(NDESourceUploader):
    name = "dash"
    __metadata__ = {
        "src_meta": {
            "sourceInfo": {
                "description": "The NICHD Data and Specimen Hub(DASH) is a centralized resource that allows researchers to share and access de-identified data from studies funded by NICHD. DASH also serves as a portal for requesting biospecimens from selected DASH studies. DASH serves as a mechanism for NICHD-funded extramural and intramural investigators to share research data from studies in accordance with NIH Data Sharing Policies. Many of the NICHD-funded research studies also collected biospecimens that are stored in the NICHD Contracted Biorepository. To provide access to these biospecimens, DASH will store and make available to other investigators the biospecimen catalog for studies that have associated research data in DASH. By supporting data and biospecimen access through DASH, NICHD aims to accelerate scientific findings and improve human health.",
                "identifier": "NICHD DASH",
                "name": "NICHD Data and Specimen Hub (DASH)",
                "schema": get_source_data(name),
                "url": "https://dash.nichd.nih.gov/",
            }
        }
    }

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        docs = standardize_funding(docs)
        docs = standardize_data(docs)
        for doc in docs:
            yield doc
