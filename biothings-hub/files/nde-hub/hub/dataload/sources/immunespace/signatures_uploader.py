from hub.dataload.nde import NDESourceUploader
from utils import iter_ndjson, nde_upload_wrapper


class ImmunespaceSignaturesUploader(NDESourceUploader):
    name = "immunespace_signatures"
    main_source = "immunespace"
    __metadata__ = {
        "merger": "merge_struct",
        "merger_kwargs": {"aslistofdict": "includedInDataCatalog", "include": ["includedInDataCatalog"]},
        "src_meta": {
            "url": "https://immunespace.org/query/results/?ordering_tab=signatures_tab",
        },
    }

    @nde_upload_wrapper
    def load_data(self, data_folder):
        yield from (doc for doc in iter_ndjson(data_folder) if doc.get("@type") == "Inference")
