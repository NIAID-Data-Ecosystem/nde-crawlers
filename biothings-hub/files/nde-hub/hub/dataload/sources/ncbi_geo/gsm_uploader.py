from utils.extract import process_descriptions
from utils.measurement_technique_helper import process_measurement_technique
from utils.pubtator import standardize_data
from utils.topic_category_helper import add_topic_category
from utils.utils import nde_upload_wrapper

from hub.dataload.nde import NDESourceSampleUploader

from .parser import parse_gsm


class GSM_Uploader(NDESourceSampleUploader):
    name = "gsm_ncbi_geo"
    main_source = "ncbi_geo"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = parse_gsm(data_folder)
        docs = standardize_data(docs)
        docs = process_measurement_technique(docs, self.main_source)
        # docs = process_descriptions(docs)
        # docs = add_topic_category(docs, self.name)
        for doc in docs:
            yield doc
