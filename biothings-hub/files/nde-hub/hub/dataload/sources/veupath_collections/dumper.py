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


class VEuPathCollections_Dumper(DockerContainerDumper):
    SRC_NAME = "veupath_collections"
    SRC_ROOT_FOLDER = os.path.join(DATA_ARCHIVE_ROOT, SRC_NAME)
    SCHEDULE = "0 16 * * 6"  # Every Saturday at 4:00 PM
    UNCOMPRESS = True

    SRC_URLS = [
        f'docker://su07?image=nde-crawlers-{SRC_NAME}-crawler&tag=latest&path=/data/{SRC_NAME}_crawled/data.ndjson&dump_command="/home/biothings/run-api-crawler.sh"&container_name={SRC_NAME}_dumper&keep_container=false'
        # &custom_cmd="/usr/bin/wget https://s3.pgkb.org/data/annotations.zip -O /tmp/annotations.zip"'
    ]
