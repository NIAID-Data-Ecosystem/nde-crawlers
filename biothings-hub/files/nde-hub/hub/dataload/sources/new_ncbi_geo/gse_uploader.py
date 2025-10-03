from hub.dataload.nde import NDESourceSampleUploader

from .parser import parse_gse


class GSE_Uploader(NDESourceSampleUploader):
    name = "gse_ncbi_geo"
    main_source = "new_ncbi_geo"

    def load_data(self, data_folder):
        docs = parse_gse(data_folder)
        for doc in docs:
            yield doc

