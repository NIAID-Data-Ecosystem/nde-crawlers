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


class DataverseDumper(DockerContainerDumper):
    SRC_NAME = "dataverse"
    SRC_ROOT_FOLDER = os.path.join(DATA_ARCHIVE_ROOT, SRC_NAME)
    SCHEDULE = "0 17 1 1,4,7,10 *"  # 1st day of Jan, Apr, Jul, Oct at 17:00 UTC
    UNCOMPRESS = True
    # NAMED_VOLUMES = [{"name": "dataverse_data", "driver":"local", "driver_opts": {"type": "none", "o": "bind", "device": "/opt/nde/data_test_1"}},
    #                  {"name": "dataverse_cache", "driver":"local", "driver_opts": {"type": "none", "o": "bind", "device": "/opt/nde/cache_test"}}]
    # VOLUMES = {"dataverse_data":{"bind": "/data"}, "dataverse_cache": {"bind": "/cache"}}
    VOLUMES = {"/opt/nde/cache": {"bind": "/cache", "mode": "rw"}}
    SRC_URLS = [
        # f'docker://su07?image=nde-crawlers_{SRC_NAME}-crawler&tag=latest&path=/data/{SRC_NAME}_crawled/data.ndjson&dump_command="/home/biothings/run-api-crawler.sh"&container_name={SRC_NAME}_dumper&volumes={json.dumps(VOLUMES)}&named_volumes={json.dumps(NAMED_VOLUMES)}'
        f'docker://su07?image=nde-crawlers_{SRC_NAME}-crawler&tag=latest&path=/data/{SRC_NAME}_crawled/data.ndjson&dump_command="/home/biothings/run-api-crawler.sh"&container_name={SRC_NAME}_dumper&volumes={json.dumps(VOLUMES)}'
        # &custom_cmd="/usr/bin/wget https://s3.pgkb.org/data/annotations.zip -O /tmp/annotations.zip"'
    ]
