from hub.dataload.nde import NDESourceSampleUploader
from utils.measurement_technique_helper import process_measurement_technique

from .parser import parse_sex


class BiosampleUploader(NDESourceSampleUploader):
    name = "biosample"

    def load_data(self, data_folder):
        docs = parse_sex(data_folder)
        docs = process_measurement_technique(docs, self.name)
        for doc in docs:
            yield doc
