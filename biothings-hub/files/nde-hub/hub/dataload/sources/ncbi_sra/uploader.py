from hub.dataload.nde import NDESourceUploader
from utils.csv_helper import get_source_data


class NCBI_SRA_Uploader(NDESourceUploader):
    name = "ncbi_sra"
    __metadata__ = {
        "src_meta": {
            "sourceInfo": {
                "description": "Sequence Read Archive(SRA) data, available through multiple cloud providers and NCBI servers, is the largest publicly available repository of high throughput sequencing data. The archive accepts data from all branches of life as well as metagenomic and environmental surveys. SRA stores raw sequencing data and alignment information to enhance reproducibility and facilitate new discoveries through data analysis.",
                "identifier": "NCBI SRA",
                "name": "NCBI SRA",
                "schema": get_source_data(name),
                "url": "https://www.ncbi.nlm.nih.gov/sra",
            }
        }
    }
