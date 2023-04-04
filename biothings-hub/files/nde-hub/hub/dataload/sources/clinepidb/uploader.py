from hub.dataload.nde import NDESourceUploader
from utils.pmid_helper import load_pmid_ctfd
from utils.csv_helper import get_source_data
from utils.utils import check_schema

class ClinEpiDB_Uploader(NDESourceUploader):
    name = "clinepidb"
    __metadata__ = {
        "src_meta": {
            "sourceInfo": {
                "description": "ClinEpiDB, launched in February 2018, is an open-access exploratory data analysis platform. We integrate data from high quality epidemiological studies, and offer tools and visualizations to explore the data within the browser in a point and click interface. We enable investigators to maximize the utility and reach of their data and to make optimal use of data released by others. ClinEpiDB is led by a team of scientists and developers based at the University of Pennsylvania, the University of Georgia, Imperial College London, and several other academic institutions. Currently, we are funded by the Bill and Melinda Gates Foundation for resource development and data integration, and by NIAID for integration of data from the International Centers of Excellence in Malaria Research (ICEMR).",
                "identifier": "ClinEpiDB",
                "name": "ClinEpiDB",
                "schema": get_source_data(name),
                "url": "https://clinepidb.org/ce/app"
            }
        }
    }

    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        for doc in docs:
            # check schema
            check_schema(doc)
            yield doc
