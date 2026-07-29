from hub.dataload.nde import NDESourceSampleUploader
from utils import nde_upload_wrapper

from .parser import parse_gsm


class GSM_Uploader(NDESourceSampleUploader):
    name = "gsm_ncbi_geo"
    main_source = "ncbi_geo"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        yield from parse_gsm(data_folder)
