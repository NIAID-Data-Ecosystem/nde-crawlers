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


class AccessClinicalDataDumper(DockerContainerDumper):
    SRC_NAME = "acd_niaid"
    SRC_ROOT_FOLDER = os.path.join(DATA_ARCHIVE_ROOT, SRC_NAME)
    SCHEDULE = None
    UNCOMPRESS = True
    NAMED_VOLUMES = {"name": "myvolume", "driver":"local", "driver_opts": {"type": "none", "o": "bind", "device": "/opt/nde/data_test_1"}}
    NAMED_VOLUMES = json.dumps(NAMED_VOLUMES)
    VOLUMES = {"myvolume":{"bind": "/data"}}
    VOLUMES = json.dumps(VOLUMES)
    SRC_URLS = [
        f'docker://su07?image=nde-crawlers_{SRC_NAME}-crawler&tag=latest&path=/data/{SRC_NAME}_crawled/data.ndjson&dump_command="/home/biothings/run-api-crawler.sh"&container_name={SRC_NAME}_dumper&keep_container=false&volumes={VOLUMES}&named_volumes={NAMED_VOLUMES}'
        # &custom_cmd="/usr/bin/wget https://s3.pgkb.org/data/annotations.zip -O /tmp/annotations.zip"'
    ]
