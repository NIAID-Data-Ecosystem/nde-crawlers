from hub.dataload.nde import NDESourceUploader
from utils.extract import process_descriptions
from utils.lineage import process_lineage
from utils.pubtator import standardize_data
from utils.topic_category_helper import add_topic_category
from utils.utils import nde_upload_wrapper


class NCBI_SRA_Uploader(NDESourceUploader):
    name = "ncbi_sra"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = standardize_data(data_folder)
        docs = process_descriptions(docs)
        docs = process_lineage(docs)
        docs = add_topic_category(docs, self.name)
        for doc in docs:
            yield doc
