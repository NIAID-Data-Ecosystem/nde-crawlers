# We can remove this line below after switched to biothings 0.11.x branch
# and remove the config_hub.py file completely
from config_hub import *
from biothings.utils.loggers import setup_default_log
import os
import urllib.parse

mongo_host = os.environ.get('MONGO_HOST', 'localhost:27017')
mongo_parsed = urllib.parse.urlparse(f'mongodb://{mongo_host}')
data_folder = os.environ.get('DATA_FOLDER', '/data/nde-hub')

# create directories in case volume changed
for _ in ['datasources', 'plugins', 'dataupload', 'diff', 'logs', 'release',
          'cache', 'run', 'esbackup']:
    os.makedirs(f'{data_folder}/{_}', exist_ok=True)

# hijack the loading progress so that the proper keys are in place before running hub
ssh_key_file = f'{data_folder}/ssh_host_key'
if not os.path.isfile(ssh_key_file):
    if os.path.exists(ssh_key_file):
        raise FileExistsError(
            f"{ssh_key_file} exists but is not a regular file")

    from cryptography.hazmat.primitives.asymmetric import rsa as pk
    from cryptography.hazmat.primitives import serialization as crypto_ser

    print("Generating SSH Keys for BioThings Hub...")
    privkey = pk.generate_private_key(65537, 2048)
    with open(ssh_key_file, 'wb') as f:
        f.write(
            privkey.private_bytes(
                crypto_ser.Encoding.PEM,
                crypto_ser.PrivateFormat.OpenSSH,
                crypto_ser.NoEncryption()
            )
        )
    pubkey = privkey.public_key().public_bytes(crypto_ser.Encoding.OpenSSH,
                                               crypto_ser.PublicFormat.OpenSSH)
    with open(f'{ssh_key_file}.pub', 'wb') as f:
        f.write(pubkey)
    print("SSH Key has been generated, Public Key:\n")
    print(pubkey.decode('ASCII'))
    print()
    del privkey, pubkey, crypto_ser, pk

DATA_ARCHIVE_ROOT = f'{data_folder}/datasources'
DATA_PLUGIN_FOLDER = f'{data_folder}/plugins'
DATA_UPLOAD_FOLDER = f'{data_folder}/dataupload'

DIFF_PATH = f"{data_folder}/diff"
RELEASE_PATH = f"{data_folder}/release"
CACHE_FOLDER = f"{data_folder}/cache"
ES_BACKUPS_FOLDER = f"{data_folder}/esbackup"

LOG_FOLDER = f"{data_folder}/logs"
logger = setup_default_log("hub", LOG_FOLDER)

RUN_DIR = f'{data_folder}/run'

DATA_SRC_SERVER = mongo_parsed.hostname or 'localhost'
DATA_SRC_PORT = mongo_parsed.port or 27017
DATA_SRC_DATABASE = 'nde_hub_src'
DATA_SRC_SERVER_USERNAME = urllib.parse.unquote(mongo_parsed.username) \
    if mongo_parsed.username else ''
DATA_SRC_SERVER_PASSWORD = urllib.parse.unquote(mongo_parsed.password) \
    if mongo_parsed.password else ''

DATA_TARGET_SERVER = DATA_SRC_SERVER
DATA_TARGET_PORT = DATA_SRC_PORT
DATA_TARGET_DATABASE = 'nde_hub'
DATA_TARGET_SERVER_USERNAME = DATA_SRC_SERVER_USERNAME
DATA_TARGET_SERVER_PASSWORD = DATA_SRC_SERVER_PASSWORD

# FIXME: deal with other uri later
DATA_HUB_DB_DATABASE = 'nde_hub_hubdb'
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
HUB_SSH_PORT = 19522
HUB_API_PORT = 19580
READONLY_HUB_API_PORT = 19581

# Hub name/icon url/version, for display purpose
HUB_NAME = "NDE Data Hub"
HUB_ICON = "https://biothings.io/static/img/sdk-icon.svg"
HUB_VERSION = "master"

USE_RELOADER = True  # so no need to restart hub when a datasource has changed

MAX_QUEUED_JOBS = 1

# INDEX_CONFIG = {
#     "indexer_select": {
#         # default
#     },
#     "env": {
#         "localhub": {
#             "host": os.environ.get('ES_HOST', 'localhost:9200'),
#             "indexer": {
#                 "args": {
#                     "timeout": 300,
#                     "retry_on_timeout": True,
#                     "max_retries": 10,
#                 },
#             },
#         },
#     },
# }

# cleanup config namespace
del os, urllib, setup_default_log, mongo_host, mongo_parsed

# This is the root data folder from output files NDEFileSystemDumper
# is monitoring, where all crawlers save output files
CRAWLER_OUTPUT_DATA_ROOT = '/data'

# Allow to override any default settings with a config_local.py file
try:
    from config_local import *
except ImportError:
    pass
