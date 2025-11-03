from hub.dataload.nde import NDESourceUploader
from utils.corrections import corrections
from utils.measurement_technique_helper import process_measurement_technique
from utils.pubtator import standardize_data
from utils.topic_category_helper import add_topic_category
from utils.utils import nde_upload_wrapper


class ARK_Uploader(NDESourceUploader):
    name = "ark"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = standardize_data(data_folder)
        docs = process_measurement_technique(docs, self.name)
        docs = corrections(docs)
        docs = add_topic_category(docs, self.name)
        for doc in docs:
            yield doc
