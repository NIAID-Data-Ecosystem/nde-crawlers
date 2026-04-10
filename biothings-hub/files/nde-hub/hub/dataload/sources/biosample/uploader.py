from hub.dataload.nde import NDESourceSampleUploader

from .parser import parse_sex


class BiosampleUploader(NDESourceSampleUploader):
    name = "biosample"

    def load_data(self, data_folder):
        docs = parse_sex(data_folder)
        for doc in docs:
            yield doc
