import os.path

import biothings
import config
from biothings.hub.dataload.dumper import DockerContainerDumper
from urllib3 import disable_warnings
from urllib3.exceptions import InsecureRequestWarning

biothings.config_for_app(config)
from config import DATA_ARCHIVE_ROOT

# Disable InsecureRequestWarning: Unverified HTTPS request is being made to host
disable_warnings(InsecureRequestWarning)


class DryadDockerDumper(DockerContainerDumper):
    SRC_NAME = "dryad_docker"
    SRC_ROOT_FOLDER = os.path.join(DATA_ARCHIVE_ROOT, SRC_NAME)
    SCHEDULE = None
    UNCOMPRESS = True
    SRC_URLS = [
        'docker://localhost?image=dryad_dumper&tag=latest&path=/data/dryad_crawled/data.ndjson&dump_command="/home/biothings/run-spider.sh"&container_name=dyrad_dumper&keep_container=false'
    ]
    __metadata__ = {"src_meta": {}}
