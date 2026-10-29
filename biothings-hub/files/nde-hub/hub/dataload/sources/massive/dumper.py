import json
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


class MassiveDumper(DockerContainerDumper):
    SRC_NAME = "massive"
    SRC_ROOT_FOLDER = os.path.join(DATA_ARCHIVE_ROOT, SRC_NAME)
    SCHEDULE = "0 16 * * 6"  # Every Saturday at 4:00 PM
    UNCOMPRESS = True
    # Keep the crawler's SQL and dataset-page caches on the host. In the container's own
    # volume they are deleted with the container after every successful dump.
    VOLUMES = {"/opt/nde/cache": {"bind": "/cache", "mode": "rw"}}
    SRC_URLS = [
        f'docker://su07?image=nde-crawlers-{SRC_NAME}-crawler&tag=latest&path=/data/{SRC_NAME}_crawled/data.ndjson&dump_command="/home/biothings/run-api-crawler.sh"&container_name={SRC_NAME}_dumper&keep_container=false&volumes={json.dumps(VOLUMES)}'
        # &custom_cmd="/usr/bin/wget https://s3.pgkb.org/data/annotations.zip -O /tmp/annotations.zip"'
    ]
