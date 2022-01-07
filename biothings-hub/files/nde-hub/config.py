from config_hub import *
from biothings.utils.loggers import setup_default_log
import os
import urllib.parse

mongo_host = os.environ.get('MONGO_HOST', 'mongodb:27017')
mongo_parsed = urllib.parse.urlparse(f'mongodb://{mongo_host}')

# create directories in case volume changed
for _ in ['datasources', 'plugins', 'dataupload', 'diff', 'logs', 'release',
          'cache', 'run', 'esbackup']:
    os.makedirs(f'/data/nde-hub/{_}', exist_ok=True)

# hijack the loading progress so that the proper keys are in place before running hub
if not os.path.isfile('/data/nde-hub/ssh_host_key'):
    if os.path.exists('/data/nde-hub/ssh_host_key'):
        raise FileExistsError(
            "/data/nde-hub/ssh_host_key exists but is not a regular file")

    from cryptography.hazmat.primitives.asymmetric import rsa as pk
    from cryptography.hazmat.primitives import serialization as crypto_ser

    print("Generating SSH Keys for BioThings Hub...")
    privkey = pk.generate_private_key(65537, 2048)
    with open('/data/nde-hub/ssh_host_key', 'wb') as f:
        f.write(
            privkey.private_bytes(
                crypto_ser.Encoding.PEM,
                crypto_ser.PrivateFormat.OpenSSH,
                crypto_ser.NoEncryption()
            )
        )
    pubkey = privkey.public_key().public_bytes(crypto_ser.Encoding.OpenSSH,
                                               crypto_ser.PublicFormat.OpenSSH)
    with open('/data/nde-hub/ssh_host_key.pub', 'wb') as f:
        f.write(pubkey)
    print("SSH Key has been generated, Public Key:\n")
    print(pubkey.decode('ASCII'))
    print()
    del privkey, pubkey, crypto_ser, pk

DATA_ARCHIVE_ROOT = f'/data/nde-hub/datasources'
DATA_PLUGIN_FOLDER = f'/data/nde-hub/plugins'
DATA_UPLOAD_FOLDER = f'/data/nde-hub/dataupload'

DIFF_PATH = f"/data/nde-hub/diff"
RELEASE_PATH = f"/data/nde-hub/release"
CACHE_FOLDER = f"/data/nde-hub/cache"
ES_BACKUPS_FOLDER = f"/data/nde-hub/esbackup"

LOG_FOLDER = f"/data/nde-hub/logs"
logger = setup_default_log("hub", LOG_FOLDER)

RUN_DIR = f'/data/nde-hub/run'

DATA_SRC_SERVER = mongo_parsed.hostname or 'localhost'
DATA_SRC_PORT = mongo_parsed.port or 27017
DATA_SRC_DATABASE = 'nde_hub_sources'
DATA_SRC_SERVER_USERNAME = urllib.parse.unquote(mongo_parsed.username) \
    if mongo_parsed.username else ''
DATA_SRC_SERVER_PASSWORD = urllib.parse.unquote(mongo_parsed.password) \
    if mongo_parsed.password else ''

DATA_TARGET_SERVER = DATA_SRC_SERVER
DATA_TARGET_PORT = DATA_SRC_PORT
DATA_TARGET_DATABASE = 'nde_hub_target'
DATA_TARGET_SERVER_USERNAME = DATA_SRC_SERVER_USERNAME
DATA_TARGET_SERVER_PASSWORD = DATA_SRC_SERVER_PASSWORD

# FIXME: deal with other uri later
DATA_HUB_DB_DATABASE = 'nde_hub_db'
HUB_DB_BACKEND = {
    "module": "biothings.utils.mongo",
    "uri": mongo_parsed.geturl(),
}

CONFIG_READONLY = False

# FIXME: deal with the version issue
# At least for BIOTHINGS_VERSION, it is trying to use the git repo version
# which does not seem to make sense
BIOTHINGS_VERSION = 'master'

# SSH port for hub console
HUB_SSH_PORT = 7022
HUB_API_PORT = 7080

# Hub name/icon url/version, for display purpose
HUB_NAME = f"Studio for NDE"
HUB_ICON = f"http://biothings.io/static/img/mygene-logo-shiny.svg"
HUB_VERSION = "master"

USE_RELOADER = True  # so no need to restart hub when a datasource has changed

MAX_QUEUED_JOBS = 1

INDEX_CONFIG = {
    "indexer_select": {
        # default
    },
    "env": {
        "localhub": {
            "host": os.environ.get('ES_HOST', 'elasticsearch:9200'),
            "indexer": {
                "args": {
                    "timeout": 300,
                    "retry_on_timeout": True,
                    "max_retries": 10,
                },
            },
        },
    },
}

# cleanup config namespace
del os, urllib, setup_default_log, mongo_host, mongo_parsed
