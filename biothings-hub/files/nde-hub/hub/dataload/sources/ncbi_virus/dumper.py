import os.path

import biothings
import config
from biothings.hub.dataload.dumper import DockerContainerDumper
from urllib3 import disable_warnings
from urllib3.exceptions import InsecureRequestWarning

biothings.config_for_app(config)
from config import DATA_ARCHIVE_ROOT

disable_warnings(InsecureRequestWarning)


class NCBI_VIRUS_Dumper(DockerContainerDumper):
    SRC_NAME = "ncbi_virus"
    SRC_ROOT_FOLDER = os.path.join(DATA_ARCHIVE_ROOT, SRC_NAME)
    SCHEDULE = "0 20 * * 6"  # Every Saturday at 8:00 PM
    UNCOMPRESS = True
    _docker_url = (
        f"docker://su07?image=nde-crawlers-{SRC_NAME}-crawler&tag=latest"
        f"&path=/data/{SRC_NAME}_crawled/data.ndjson"
        '&dump_command="/home/biothings/run-api-crawler.sh"'
        f"&container_name={SRC_NAME}_dumper&keep_container=false"
    )
    SRC_URLS = [
        _docker_url,
    ]
