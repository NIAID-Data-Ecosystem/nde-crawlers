from hub.dataload.nde import NDESourceUploader
from utils.csv_helper import get_source_data


class LINCSUploader(NDESourceUploader):
    name = "lincs"
    __metadata__ = {
        "src_meta": {
            "sourceInfo": {
                "description": "The BD2K-LINCS DCIC is comprised of four major components: Integrated Knowledge Environment (IKE), Data Science Research (DSR), Community Training and Outreach (CTO) and Consortium Coordination and Administration (CCA). The Center is constructing a high-capacity scalable integrated knowledge environment enabling federated access, intuitive querying and integrative analysis and visualization across all LINCS resources and many additional external data types from other relevant resources. The Center’s data science research projects are aimed at addressing various data integration and intracellular molecular regulatory network challenges. The Center aims to develop: 1) methods to connect cellular and organismal phenotypes with molecular cellular signatures, and 2) novel data visualization methods for dynamically interacting with large-genomics and proteomics datasets.",
                "identifier": "LINCS",
                "name": "BD2K-LINCS DCIC",
                "schema": get_source_data(name),
                "url": "https://lincsportal.ccs.miami.edu/"
            }
        }
    }
