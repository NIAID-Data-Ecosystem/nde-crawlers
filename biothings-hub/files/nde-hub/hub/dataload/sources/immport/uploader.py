from hub.dataload.nde import NDESourceUploader
from utils.corrections import corrections
from utils.disambiguating_description import add_disambiguating_description
from utils.extract import process_descriptions
from utils.funding_helper import standardize_funding
from utils.pubtator import standardize_data
from utils.topic_category_helper import add_topic_category
from utils.utils import nde_upload_wrapper


class ImmPortUploader(NDESourceUploader):
    name = "immport"
    __metadata__ = {
        "merger": "merge_struct",
        "merger_kwargs": {"aslistofdict": "includedInDataCatalog", "include": ["includedInDataCatalog"]},
        "src_meta": {
            "url": "https://www.immport.org/shared/home",
            "license_url": "https://docs.immport.org/home/agreement/",
        },
    }

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = standardize_data(data_folder)
        docs = standardize_funding(docs)
        docs = process_descriptions(docs)
        docs = corrections(docs)
        docs = add_disambiguating_description(docs, self.name)
        docs = add_topic_category(docs, self.name)
        for doc in docs:
            yield doc
