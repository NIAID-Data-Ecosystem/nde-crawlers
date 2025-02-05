from hub.dataload.nde import NDESourceUploader
from utils.corrections import corrections
from utils.extract import process_descriptions
from utils.funding_helper import standardize_funding
from utils.pmid_helper import load_pmid_ctfd
from utils.pubtator import standardize_data
from utils.topic_category_helper import add_topic_category
from utils.utils import nde_upload_wrapper


class OmicsDIUploader(NDESourceUploader):
    name = "omicsdi"
    __metadata__ = {
        "src_meta": {"url": "https://www.omicsdi.org/search", "license_url": "https://www.ebi.ac.uk/licencing"},
        "merger": "merge_struct",
        "merger_kwargs": {"aslistofdict": "includedInDataCatalog", "include": ["includedInDataCatalog"]},
    }

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        docs = standardize_funding(docs)
        docs = standardize_data(docs)
        docs = process_descriptions(docs)
        docs = corrections(docs)
        docs = add_topic_category(docs, self.name)
        for doc in docs:
            yield doc
