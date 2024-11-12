from hub.dataload.nde import NDESourceUploader
from utils.extract import process_descriptions
from utils.funding_helper import standardize_funding
from utils.measurement_technique_helper import process_measurement_technique
from utils.pmid_helper import load_pmid_ctfd
from utils.pubtator import standardize_data
from utils.topic_category_helper import add_topic_category
from utils.utils import nde_upload_wrapper


class ReframedbUploader(NDESourceUploader):
    name = "reframedb"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        docs = standardize_funding(docs)
        docs = standardize_data(docs)
        docs = process_measurement_technique(docs, self.name)
        docs = process_descriptions(docs)
        docs = add_topic_category(docs, self.name)
        for doc in docs:
            yield doc
