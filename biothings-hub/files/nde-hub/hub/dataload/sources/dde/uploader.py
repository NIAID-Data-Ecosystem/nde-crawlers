from hub.dataload.nde import NDESourceUploader
from utils.corrections import corrections
from utils.extract import process_descriptions
from utils.funding_helper import standardize_funding
from utils.in_defined_term_set import handle_dde_docs
from utils.lineage import process_lineage
from utils.pmid_helper import load_pmid_ctfd, standardize_fields
from utils.pubtator import standardize_data
from utils.topic_category_helper import add_topic_category
from utils.utils import nde_upload_wrapper


class DDEUploader(NDESourceUploader):
    name = "dde"
    __metadata__ = {
        "src_meta": {
            "url": "https://discovery.biothings.io/api/dataset/",
            "license_url": "https://creativecommons.org/licenses/by/4.0/",
            "license": "Creative Commons Attribution 4.0 International",
        }
    }

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        docs = standardize_fields(docs)
        docs = handle_dde_docs(docs)
        docs = standardize_funding(docs)
        docs = standardize_data(docs)
        docs = process_descriptions(docs)
        docs = process_lineage(docs)
        docs = corrections(docs)
        docs = add_topic_category(docs, self.name)
        for doc in docs:
            yield doc
