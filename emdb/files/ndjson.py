import datetime
import logging
import os
import platform
import traceback

import orjson

import emdb

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

# set the release string to be ISO date format
# minute precision is good enough but feel free to change
release_string = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
dirname = os.path.join("/data", "emdb_crawled")
os.makedirs(dirname, exist_ok=True)
release_filename = os.path.join(dirname, "release.txt")
# so that updates are as atomic as possible, using rename
final_data_filename = os.path.join(dirname, "data.ndjson")
tmp_filename = f"{final_data_filename}.{platform.node()}.{os.getpid()}"
rl_tmp_filename = f"{release_filename}.{platform.node()}.{os.getpid()}"
with open(rl_tmp_filename, "w") as release_file:
    release_file.write(release_string)
os.makedirs(os.path.dirname(tmp_filename), exist_ok=True)
fd = open(tmp_filename, "wb")

is_parsed = False
# run parser
try:
    docs = emdb.parse()
    for doc in docs:
        line = orjson.dumps(doc) + b"\n"
        fd.write(line)
    is_parsed = True
# parser failed
except Exception as e:
    fd.close()
    logger.warning("Errors occurred while running, so not saving potentially corrupt data.")
    os.unlink(tmp_filename)
    os.unlink(rl_tmp_filename)
    logger.error(e)
    logger.error(traceback.format_exc())
finally:
    fd.close()


if is_parsed:
    try:
        # it is critical that we use rename here so that
        # the files have a different inode
        # because NDE Dumper uses hardlinks
        os.rename(tmp_filename, final_data_filename)
        os.rename(rl_tmp_filename, release_filename)
    except Exception:
        # cleanup the mess by attempting to remove everything
        # so that bad data does not propagate
        logger.error("Error updating data/release, will remove corrupt files")
        errors = []
        for fn in [tmp_filename, final_data_filename, rl_tmp_filename, release_filename]:
            try:
                os.unlink(fn)
            except FileNotFoundError:
                pass
            except Exception as e:
                logger.critical(f"Failed to remove corrupt file {fn}")
                errors.append(e)
        if errors:
            raise RuntimeError(errors)
