"""PubTator3 species and disease annotation dumps for the citations stage.

https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/README.txt
"""

import os

import biothings
import config

biothings.config_for_app(config)
from biothings.hub.dataload.dumper import LastModifiedHTTPDumper
from config import DATA_ARCHIVE_ROOT


class PubTatorDumper(LastModifiedHTTPDumper):
    SRC_NAME = "pubtator"
    SRC_ROOT_FOLDER = os.path.join(DATA_ARCHIVE_ROOT, SRC_NAME)
    SRC_URLS = [
        "https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/species2pubtator3.gz",
        "https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/disease2pubtator3.gz",
    ]
    SCHEDULE = "0 6 * * 1"  # Mondays at 06:00; PubTator3 updates monthly
    ARCHIVE = False  # both files stay in one folder when only one changes
