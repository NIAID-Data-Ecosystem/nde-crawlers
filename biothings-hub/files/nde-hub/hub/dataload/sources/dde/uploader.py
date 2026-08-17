from hub.dataload.nde import NDESourceUploader
from utils import iter_ndjson, nde_upload_wrapper
from utils.citations import standardize_fields
from utils.dde import handle_dde_docs


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
        # DDE submitters pick their own ontology terms, so curate those before
        # the pipeline standardizes whatever is left.
        docs = standardize_fields(iter_ndjson(data_folder))
        yield from handle_dde_docs(docs)
