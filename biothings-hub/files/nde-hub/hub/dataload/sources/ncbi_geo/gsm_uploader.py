from hub.dataload.nde import NDESourceSampleUploader

from .parser import parse_gsm


class GSM_Uploader(NDESourceSampleUploader):
    name = "gsm_ncbi_geo"
    main_source = "ncbi_geo"

    def load_data(self, data_folder):
        docs = parse_gsm(data_folder)
        for doc in docs:
            yield doc

