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


class DockStoreDockerDumper(DockerContainerDumper):
    SRC_NAME = "dockstore_docker"
    SRC_ROOT_FOLDER = os.path.join(DATA_ARCHIVE_ROOT, SRC_NAME)
    SCHEDULE = None
    UNCOMPRESS = True
    SRC_URLS = [
        'docker://localhost?image=dockstore_dumper&tag=latest&path=/data/dockstore_crawled/data.ndjson&dump_command="/home/biothings/run-dockstore.sh"&container_name=dockstore_dumper&keep_container=false'
        # &custom_cmd="/usr/bin/wget https://s3.pgkb.org/data/annotations.zip -O /tmp/annotations.zip"'
    ]
    __metadata__ = {"src_meta": {}}
