from hub.dataload.nde import NDESourceUploader

from .parser import parse_gse


class GSE_Uploader(NDESourceUploader):
    name = "gse_ncbi_geo"
    main_source = "new_ncbi_geo"

    def load_data(self, data_folder):
        docs = parse_gse(data_folder)
        for doc in docs:
            yield doc

