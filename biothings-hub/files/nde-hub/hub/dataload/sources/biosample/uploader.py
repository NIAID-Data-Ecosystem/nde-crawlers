from hub.dataload.nde import NDESourceSampleUploader
from utils.measurement_technique_helper import process_measurement_technique
from utils.pubtator import standardize_data
from utils.utils import nde_upload_wrapper

from .parser import parse_sex_docs


class BiosampleUploader(NDESourceSampleUploader):
    name = "biosample"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = standardize_data(data_folder)
        docs = parse_sex_docs(docs)
        docs = process_measurement_technique(docs, self.name)
        for doc in docs:
            yield doc
