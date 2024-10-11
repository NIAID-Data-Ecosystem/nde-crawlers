from hub.dataload.nde import NDESourceUploader
from utils.csv_helper import get_source_data
from utils.disambiguating_description import add_disambiguating_description
from utils.extract import process_descriptions
from utils.funding_helper import standardize_funding
from utils.pmid_helper import load_pmid_ctfd
from utils.pubtator import standardize_data
from utils.topic_category_helper import add_topic_category
from utils.utils import nde_upload_wrapper


class ClinEpiDB_Uploader(NDESourceUploader):
    name = "clinepidb"
    __metadata__ = {
        "src_meta": {
            "sourceInfo": {
                "description": "ClinEpiDB, launched in February 2018, is an open-access exploratory data analysis platform. We integrate data from high quality epidemiological studies, and offer tools and visualizations to explore the data within the browser in a point and click interface. We enable investigators to maximize the utility and reach of their data and to make optimal use of data released by others. ClinEpiDB is led by a team of scientists and developers based at the University of Pennsylvania, the University of Georgia, Imperial College London, and several other academic institutions. Currently, we are funded by the Bill and Melinda Gates Foundation for resource development and data integration, and by NIAID for integration of data from the International Centers of Excellence in Malaria Research (ICEMR).",
                "identifier": "ClinEpiDB",
                "name": "ClinEpiDB",
                "schema": get_source_data(name),
                "url": "https://clinepidb.org/ce/app",
            }
        }
    }

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        docs = standardize_funding(docs)
        docs = standardize_data(docs)
        docs = process_descriptions(docs)
        docs = add_topic_category(docs, self.name)
        docs = add_disambiguating_description(docs, self.name)
        for doc in docs:
            yield doc
