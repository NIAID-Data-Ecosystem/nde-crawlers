from hub.dataload.nde import NDESourceUploader
from utils import nde_upload_wrapper

from .parser import parse_gse


class GSE_Uploader(NDESourceUploader):
    name = "gse_ncbi_geo"
    main_source = "ncbi_geo"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        yield from parse_gse(data_folder)
