import datetime
import logging
import os
import platform
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Manager, Process

import orjson

from omicsdi import get_dataset_names, process_dataset_records

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("nde-logger")


def write_ndjson(queue):
    """
    Writes documents from a generator to an NDJSON file.
    """
    # Set up file paths
    release_string = datetime.datetime.now(datetime.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ"
    )
    dirname = os.path.join("/data", "omicsdi_crawled")
    os.makedirs(dirname, exist_ok=True)
    release_filename = os.path.join(dirname, "release.txt")
    final_data_filename = os.path.join(dirname, "data.ndjson")
    tmp_filename = f"{final_data_filename}.{platform.node()}.{os.getpid()}"
    rl_tmp_filename = f"{release_filename}.{platform.node()}.{os.getpid()}"

    # Write release string to a temporary file
    with open(rl_tmp_filename, "w") as release_file:
        release_file.write(release_string)
    os.makedirs(os.path.dirname(tmp_filename), exist_ok=True)
    fd = open(tmp_filename, "wb")

    is_parsed = False
    try:
        while True:
            record = queue.get()
            if record == "STOP":
                break
            if record:
                line = orjson.dumps(record) + b"\n"
                fd.write(line)
        is_parsed = True
    except Exception:
        logger.error(
            "Errors occurred while running, so not saving potentially corrupt data."
        )
        logger.error(traceback.format_exc())
        fd.close()
        os.unlink(tmp_filename)
        os.unlink(rl_tmp_filename)
        return

    fd.close()

    if is_parsed:
        try:
            os.rename(tmp_filename, final_data_filename)
            os.rename(rl_tmp_filename, release_filename)
        except Exception:
            logger.error("Error updating data/release, will remove corrupt files")
            errors = []
            for fn in [
                tmp_filename,
                final_data_filename,
                rl_tmp_filename,
                release_filename,
            ]:
                try:
                    os.unlink(fn)
                except FileNotFoundError:
                    pass
                except Exception as e:
                    logger.critical(f"Failed to remove corrupt file {fn}")
                    errors.append(e)
            if errors:
                raise RuntimeError(errors)


if __name__ == "__main__":

    dataset_names = get_dataset_names()

    # Use Manager to create a queue that can be shared across processes
    with Manager() as manager:
        queue = manager.Queue()
        # Start the writer process
        writer_process = Process(target=write_ndjson, args=(queue,))
        writer_process.start()
        logger.info("Writer process started.")

        # Use ProcessPoolExecutor with a maximum of 4 workers
        with ProcessPoolExecutor(max_workers=4) as executor:
            future_to_dataset = {
                executor.submit(process_dataset_records, name, queue): name
                for name in dataset_names
            }

            # Track the progress of each process
            for future in as_completed(future_to_dataset):
                dataset_name = future_to_dataset[future]
                try:
                    future.result()  # Wait for the process to complete
                    logger.info(f"Process completed for dataset: {dataset_name}")
                except Exception as e:
                    logger.error(
                        f"Process failed for dataset: {dataset_name} with error: {e}"
                    )

        # Signal the generator to stop
        queue.put("STOP")
