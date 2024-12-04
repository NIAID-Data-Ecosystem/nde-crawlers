import ftplib
import functools
import json
import logging
import os
import re
import time
import traceback

import xmltodict
from lxml import etree, html

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("nde-logger")

# Define the FTP URL
ftp_url = "ftp.ncbi.nlm.nih.gov"
root_directory = "/dbgap/studies"


def retry(retry_num, retry_sleep_sec):
    """
    retry help decorator.
    :param retry_num: the retry num; retry sleep sec
    :return: decorator
    """

    def decorator(func):
        """decorator"""

        # preserve information about the original function, or the func name will be "wrapper" not "func"
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            """wrapper"""
            for attempt in range(retry_num):
                try:
                    return func(*args, **kwargs)  # should return the raw function's return value
                except (BrokenPipeError, EOFError) as err:
                    logger.error(err)
                    logger.error(traceback.format_exc())
                    logger.error(f"Broken pipe encountered. Retrying in {retry_sleep_sec} seconds...")
                    time.sleep(retry_sleep_sec)
                    ftp = args[0]
                    ftp.connect("ftp.ncbi.nlm.nih.gov")  # Reconnect to the FTP server
                    ftp.login()
                except Exception as err:
                    logger.error(err)
                    logger.error(traceback.format_exc())
                    logger.error(f"Rate limit encountered. Retrying in {retry_sleep_sec} seconds...")
                    time.sleep(retry_sleep_sec)

                logger.info(
                    "Retrying failed func %s. Trying attempt %s of %s.",
                    func.__name__,
                    attempt + 1,
                    retry_num,
                )
            logger.error("func %s retry failed", func.__name__)
            raise Exception("Exceed max retry num: {} failed".format(retry_num))

        return wrapper

    return decorator


def get_highest_version(subdirs):
    highest_version = (0, 0)
    highest_dir = ""
    for subdir in subdirs:
        match = re.match(r"phs\d+\.v(\d+)\.p(\d+)", subdir)
        if match:
            version = int(match.group(1))
            part = int(match.group(2))
            if (version, part) > highest_version:
                highest_version = (version, part)
                highest_dir = subdir
    return highest_dir


@retry(retry_num=5, retry_sleep_sec=5)
def download_highest_version_xml(ftp, phs_dir, highest_version_dir, download_dir="./xml"):
    ftp.cwd(highest_version_dir)
    files = ftp.nlst()

    highest_version = (0, 0)
    highest_xml = ""
    for file in files:
        match = re.match(r"GapExchange_phs\d+\.v(\d+)\.p(\d+)\.xml", file)
        if match:
            version = int(match.group(1))
            part = int(match.group(2))
            if (version, part) > highest_version:
                highest_version = (version, part)
                highest_xml = file

    if highest_xml:
        if not os.path.exists(download_dir):
            os.makedirs(download_dir)
        local_filepath = os.path.join(download_dir, phs_dir + ".xml")
        with open(local_filepath, "wb") as f:
            ftp.retrbinary(f"RETR {highest_xml}", f.write)
        logger.info(f"Downloaded: {highest_xml} to {local_filepath}")
    else:
        logger.info(f"No xml file found in {highest_version_dir}")


def download():
    # Ensure the download directory exists
    download_dir = "./xml"
    os.makedirs(download_dir, exist_ok=True)
    json_dir = "./json"
    os.makedirs(json_dir, exist_ok=True)

    # Use the context manager to manage the FTP connection
    with ftplib.FTP(ftp_url) as ftp:
        ftp.login()  # Anonymous login
        ftp.cwd(root_directory)

        # List all directories
        directories = ftp.nlst()

        # Filter directories that start with "phs"
        phs_directories = [d for d in directories if d.startswith("phs")]

        for phs in phs_directories:
            ftp.cwd(f"{root_directory}/{phs}")
            subdirs = ftp.nlst()
            highest_version = get_highest_version(subdirs)
            highest_version_dir = f"{root_directory}/{phs}/{highest_version}"
            if highest_version_dir:
                logger.info(f"Highest version directory in {phs}: {highest_version_dir}")
                download_highest_version_xml(ftp, phs, highest_version_dir)
                local_filepath = os.path.join(download_dir, phs + ".xml")

                if os.path.exists(local_filepath):
                    with open(local_filepath, "rb") as xml_file:
                        root = html.fromstring(xml_file.read())
                        parsed_xml = xmltodict.parse(etree.tostring(root))
                        parsed_xml["gapexchange"]["studies"]["study"].pop("analyses", None)
                        parsed_xml["gapexchange"]["studies"]["study"].pop("annotations", None)
                    with open(f"{json_dir}/{phs}.json", "w", encoding="utf-8") as json_file:
                        json.dump(parsed_xml, json_file, indent=2)
                else:
                    logger.info(f"File not found: {local_filepath}")
