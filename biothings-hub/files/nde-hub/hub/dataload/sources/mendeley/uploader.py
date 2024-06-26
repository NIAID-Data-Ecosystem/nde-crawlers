from hub.dataload.nde import NDESourceUploader
from utils.topic_category_helper import add_topic_category
from utils.utils import nde_upload_wrapper


class MendeleyUploader(NDESourceUploader):
    name = "mendeley"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = add_topic_category(data_folder, self.name)
        for doc in docs:
            yield doc
