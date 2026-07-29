from hub.dataload.nde import NDESourceUploader
from utils import iter_ndjson, nde_upload_wrapper
from utils.clinical_trials import load_ct_wrapper


class VivliUploader(NDESourceUploader):
    name = "vivli"
    __metadata__ = {
        "src_meta": {
            "url": "https://search.vivli.org/",
            "license_url": "https://vivli.org/resources/vivli-data-use-agreement/",
        }
    }

    @nde_upload_wrapper
    @load_ct_wrapper
    def load_data(self, data_folder):
        yield from iter_ndjson(data_folder)
