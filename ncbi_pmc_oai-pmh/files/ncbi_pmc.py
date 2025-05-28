import datetime
import json
import logging
import os
import shutil
import socket
import tarfile
import unicodedata
from ftplib import FTP, error_temp
from time import sleep
from urllib.error import ContentTooShortError
from xml.etree import ElementTree

import wget
from dateutil.parser import parse
from lxml import etree
from sql_database import NDEDatabase

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("nde-logger")


class NCBI_PMC(NDEDatabase):
    # override variables
    SQL_DB = "ncbi_pmc.db"
    EXPIRE = datetime.timedelta(days=90)

    def get_metadata(self, record):
        """Converts the XML record to a dictionary"""
        metadata_dict = {}
        for el in record.iter():
            key = el.tag.split("}")[-1]
            if key not in metadata_dict:
                metadata_dict[key] = [x.strip() for x in el.itertext() if x.strip() != ""]
            else:
                text_list = [x.strip() for x in el.itertext() if x.strip() != ""]
                for text in text_list:
                    metadata_dict[key].append(text)
        metadata_dict.pop("article")
        return metadata_dict

    def get_zipped_files(self, last_updated=None):
        """Used to get the list of zipped files from the FTP server"""
        try_count = 0
        while True:
            try:
                logger.info("Connecting to ftp.ncbi.nlm.nih.gov")
                ftp = FTP("ftp.ncbi.nlm.nih.gov")
                ftp.login()
                logger.info("Connected to ftp.ncbi.nlm.nih.gov")
                ftp.cwd("pub/pmc/oa_bulk/oa_comm/xml/")
                logger.info("Changed directory to pub/pmc/oa_bulk/oa_comm/xml/")
                files = ftp.nlst()
                break
            except (error_temp, BrokenPipeError, socket.timeout) as e:
                try_count += 1
                if try_count > 10:
                    raise e
                else:
                    logger.error("Error getting file list. %s. Retrying...", e)
                    logger.info("Try count: %s of 10", try_count)
                    sleep(try_count * 2)
        zipped_filenames = [x for x in files if x.endswith(".tar.gz")]

        if last_updated:
            new_files = []
            for zipped_filename in zipped_filenames:
                # get the last modified date of the file
                last_modified = ftp.sendcmd("MDTM " + zipped_filename)
                # convert to YYYY-MM-DD format
                last_modified = datetime.datetime.strptime(last_modified[4:], "%Y%m%d%H%M%S").strftime("%Y-%m-%d")
                # if the file has been updated since last_updated, add it to the list of files to download
                if last_modified > last_updated and zipped_filename.endswith(".tar.gz"):
                    new_files.append(zipped_filename)
            logger.info(f"Found {len(new_files)} new files since last update on {last_updated}")
            return new_files
        else:
            logger.info(f"Found {len(zipped_filenames)} files")
            return zipped_filenames

    def load_cache(self):
        zipped_filenames = self.get_zipped_files()

        logger.info("Starting to download files")

        count = 0
        for zipped_filename in zipped_filenames:
            retry_count = 0
            while True:
                try:
                    logger.info("Downloading %s", zipped_filename)
                    wget.download("https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/xml/" + zipped_filename)
                    break
                except ContentTooShortError as e:
                    retry_count += 1
                    if retry_count > 10:
                        raise e
                    else:
                        logger.error("Error downloading file. %s. Retrying...", e)
                        logger.info("Try count: %s of 10", retry_count)
                        sleep(retry_count * 2)

            logger.info("Downloaded %s", zipped_filename)

            unzipped_file = tarfile.open(zipped_filename)
            logger.info("Unzipped %s", zipped_filename)

            unzipped_filepaths = unzipped_file.getnames()
            logger.info("Retrieved list of unzipped files in %s", zipped_filename)

            unzipped_file.extractall()
            logger.info("Extracted %s", zipped_filename)

            unzipped_file.close()
            logger.info("Closed %s", zipped_filename)

            os.remove(zipped_filename)
            logger.info("Deleted %s", zipped_filename)

            logger.info("Starting to yield files to cache")
            for xml_filepath in unzipped_filepaths:
                count += 1
                if count % 1000 == 0:
                    logger.info("Yielded %s files to cache", count)

                record_file = open(xml_filepath)
                # logger.info('Opened %s', xml_filepath)

                try:
                    record = ElementTree.parse(record_file)
                    root = record.getroot()
                    metadata_dict = self.get_metadata(root)
                    record_xml_string = ElementTree.tostring(root, encoding="unicode")

                    doc = {"metadata": metadata_dict, "xml": record_xml_string}

                    # Supplemental Data
                    supplemental_data_arr = root.findall(".//supplementary-material")
                    if len(supplemental_data_arr) == 0:
                        logger.info(
                            f"No supplemental data found for {xml_filepath.split('/')[1].replace('.xml', '')}, skipping"
                        )
                        continue

                    yield (xml_filepath.split("/")[1].replace(".xml", ""), json.dumps(doc))
                    # logger.info('Yielded %s', xml_filepath)

                    record_file.close()
                    # logger.info('Closed %s', xml_filepath)
                except ElementTree.ParseError as e:
                    logger.error("Error parsing %s", xml_filepath)
                    logger.error(e)
                    record_file.close()
                    os.remove(xml_filepath)
                    continue

            directories = set([x.split("/")[0] for x in unzipped_filepaths])
            [shutil.rmtree(x, ignore_errors=True) for x in directories]
            logger.info("Deleted directories: %s", directories)
            break

    def parse(self, records):
        no_supplementary_material_count = 0
        for record in records:
            data = json.loads(record[1])
            try:
                root = etree.fromstring(data["xml"])
            except etree.XMLSyntaxError as e:
                logger.error("Error parsing XML for %s", record[0])
                logger.error(e)
                continue
            metadata = data["metadata"]

            output = {
                "includedInDataCatalog": {
                    "@type": "DataCatalog",
                    "name": "NCBI PMC",
                    "versionDate": datetime.date.today().strftime("%Y-%m-%d"),
                    "url": "https://www.ncbi.nlm.nih.gov/pmc/",
                },
                "@type": "Dataset",
            }

            date_published = root.find('.//pub-date[@pub-type="pmc-release"]')
            if date_published is not None:
                date_published_string = ""
                year = date_published.find("year")
                month = date_published.find("month")
                day = date_published.find("day")
                if year is not None:
                    date_published_string = year.text
                if month is not None:
                    if len(month.text) == 1:
                        date_published_string += f"-0{month.text}"
                    else:
                        date_published_string += "-" + month.text
                if day is not None:
                    if len(day.text) == 1:
                        date_published_string += f"-0{day.text}"
                    else:
                        date_published_string += "-" + day.text
                if date_published_string != "":
                    try:
                        date_published_string = parse(date_published_string)
                        date_published_string = date_published_string.date()
                        date_published_string = date_published_string.isoformat()
                        output["datePublished"] = date_published_string
                    except ValueError:
                        logger.info(f"Unable to parse date: {date_published_string}")

            citation_dict = {}
            if identifiers := metadata.get("article-id"):
                for identifier in identifiers:
                    if identifier.startswith("PMC"):
                        output["identifier"] = identifier
                        output["url"] = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{identifier}"
                        output["includedInDataCatalog"][
                            "archivedAt"
                        ] = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{identifier}"
                        output["_id"] = "NCBI_PMC_" + identifier
                    elif identifier.startswith("10."):
                        citation_dict["doi"] = identifier
            if "identifier" not in output:
                logger.info("Article has been deleted, skipping")
                continue

            notes = root.find(".//notes")
            if notes is not None:
                volume = notes.find(".//volume")
                if volume is not None:
                    citation_dict["volumeNumber"] = volume.text

            # Supplemental Data
            supplemental_data_arr = root.findall(".//supplementary-material")
            distribuiton_list = []
            for supplemental_data in supplemental_data_arr:
                distribution_obj = {}
                caption = supplemental_data.find("caption")
                if caption is not None:
                    file_name = caption.find("title")
                    if file_name is not None:
                        if file_name.text is not None:
                            name = unicodedata.normalize("NFKD", file_name.text)
                            distribution_obj["name"] = name
                    description = caption.find("p")
                    if description is not None:
                        distribution_obj["description"] = description.text
                media = supplemental_data.find("media")
                if media is not None:
                    file_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{output['identifier']}/bin/" + media.get(
                        "{http://www.w3.org/1999/xlink}href"
                    )
                    distribution_obj["url"] = file_url
                if bool(distribution_obj):
                    distribuiton_list.append(distribution_obj)
            if len(distribuiton_list):
                output["distribution"] = distribuiton_list
            else:
                logger.info(f'No supplemental data found for {output["identifier"]}, skipping')
                no_supplementary_material_count += 1
                continue

            if journal_title := metadata.get("journal-title"):
                citation_dict["journalName"] = journal_title[0]

            pub_date = root.find('.//pub-date[@pub-type="epub"]')
            if pub_date is not None:
                date_string = ""
                year = pub_date.find("year")
                month = pub_date.find("month")
                day = pub_date.find("day")
                if year is not None:
                    date_string = year.text
                if month is not None:
                    if len(month.text) == 1:
                        date_string += f"-0{month.text}"
                    else:
                        date_string += "-" + month.text
                if day is not None:
                    if len(day.text) == 1:
                        date_string += f"-0{day.text}"
                    else:
                        date_string += "-" + day.text
                if date_string != "":
                    try:
                        date_string = parse(date_string)
                        date_string = date_string.date()
                        date_string = date_string.isoformat()
                        citation_dict["datePublished"] = date_string
                    except ValueError:
                        logger.info(f"Unable to parse date: {date_string}")

            # if issn := metadata.get('issn'):
            #     output['issueNumber'] = issn[0]

            # if license := metadata.get('license_ref'):
            #     output['license'] = license[0]

            if keywords := metadata.get("kwd"):
                output["keywords"] = keywords

            funding_group = root.findall(".//funding-source")
            funder_list = []
            for funder in funding_group:
                if funder.text is not None:
                    if funder.text.strip() == "":
                        continue
                    funder_list.append({"funder": {"name": funder.text}})
                else:
                    name = funder.find(".//institution")
                    if name is not None and name.text is not None:
                        if name.text.strip() == "":
                            continue
                        funder_list.append({"funder": {"name": name.text}})
            if len(funder_list):
                output["funding"] = funder_list

            publisher_sec = root.find(".//publisher")
            if publisher_sec is not None:
                publisher_list = publisher_sec.getchildren()
                for publisher in publisher_list:
                    if publisher.tag == "publisher-name":
                        output["sdPublisher"] = {"name": publisher.text}

            # if title := metadata.get('article-title'):
            #     if title[0] is not None:
            #         output['name'] = 'Supplementary materials in ' + title[0]

            article_meta_sec = root.find(".//article-meta")
            if article_meta_sec is not None:
                if "volume" not in citation_dict:
                    volume = article_meta_sec.find(".//volume")
                    if volume is not None:
                        citation_dict["volumeNumber"] = volume.text

                title_group = article_meta_sec.find(".//title-group")
                if title_group is not None:
                    article_title = title_group.find(".//article-title")
                    if article_title is not None:
                        article_name = "".join(article_title.itertext())
                        if article_name != "" or article_name != " ":
                            output["name"] = f'Supplementary materials in "{article_name}"'
                            citation_dict["name"] = article_name

            pmid = root.find(".//article-id[@pub-id-type='pmid']")
            if pmid is not None:
                citation_dict["pmid"] = pmid.text
                citation_dict["identifier"] = "PMID:" + pmid.text
                citation_dict["url"] = f"https://pubmed.ncbi.nlm.nih.gov/{pmid.text}"

            contrib_list = root.findall(".//contrib")
            author_list = []
            for contrib in contrib_list:
                author_dict = {}
                author_id = contrib.find("contrib-id")
                if author_id is not None:
                    author_dict["identifier"] = author_id.text
                surname = contrib.find(".//surname")
                if surname is not None:
                    author_dict["familyName"] = surname.text
                given_name = contrib.find(".//given-names")
                if given_name is not None:
                    author_dict["givenName"] = given_name.text
                if bool(author_dict):
                    author_dict["@type"] = "Person"
                    author_list.append(author_dict)
            if len(author_list):
                output["author"] = author_list
                citation_dict["author"] = author_list

            # Data Availability
            # data_availability_sec = root.find(
            #     './/sec[@sec-type="data-availability"]')
            # if data_availability_sec is not None:
            #     data_availability_titles = data_availability_sec.findall(
            #         'p')
            #     if data_availability_titles is not None:
            #         for title in data_availability_titles:
            #             data_availability_dict = {}
            #             content_list = []
            #             data_availability_text = title.text
            #             if data_availability_text is not None:
            #                 data_availability_dict['title'] = data_availability_text

            #             data_availability_content_titles = title.findall(
            #                 './/p')
            #             for title in data_availability_content_titles:
            #                 data_availability_content_link = title.find(
            #                     './/ext-link')
            #                 content_list.append({
            #                     'link_title': title.text.replace(' (', ''), 'link': data_availability_content_link.text})
            #             data_availability_dict['content'] = content_list

            # abstract
            abstract_secs = root.findall(".//abstract")
            if len(abstract_secs) > 0:
                description_text = ""
                for abstract_sec in abstract_secs:
                    abstract_title = abstract_sec.find(".//title")
                    abstract_text = abstract_sec.find(".//p")
                    if abstract_title is not None:
                        abstract_title_text = "".join(abstract_title.itertext())
                        if abstract_title_text is not None:
                            if description_text != "":
                                description_text += f"\n{abstract_title_text}\n"
                            else:
                                description_text += abstract_title_text
                    if abstract_text is not None:
                        abstract_text_text = "".join(abstract_text.itertext())
                        if abstract_text_text is not None:
                            if description_text != "":
                                description_text += f"\n{abstract_text_text}"
                            else:
                                description_text += abstract_text_text
                if description_text != "":
                    output["description"] = description_text

            body_sec = root.find(".//body")
            if body_sec is not None:
                description_text = ""
                sec_list = body_sec.findall(".//sec")
                if len(sec_list) > 0:
                    for sec in sec_list:
                        sec_children = sec.getchildren()
                        for sec_child in sec_children:
                            if sec_child.tag == "title":
                                sec_title_text = "".join(sec_child.itertext())
                                if sec_title_text is not None:
                                    if description_text != "":
                                        description_text += f"\n{sec_title_text}\n"
                                    else:
                                        description_text += sec_title_text
                            if sec_child.tag == "p":
                                sec_text_text = "".join(sec_child.itertext())
                                if sec_text_text is not None:
                                    if description_text != "":
                                        description_text += f"\n{sec_text_text}"
                                    else:
                                        description_text += sec_text_text
                    if description_text != "":
                        output["description"] = description_text
                else:
                    body_children = body_sec.getchildren()
                    for body_child in body_children:
                        if body_child.tag == "title":
                            body_title_text = "".join(body_child.itertext())
                            if body_title_text is not None:
                                if description_text != "":
                                    description_text += f"\n{body_title_text}\n"
                                else:
                                    description_text += body_title_text
                        if body_child.tag == "p":
                            body_text_text = "".join(body_child.itertext())
                            if body_text_text is not None:
                                if description_text != "":
                                    description_text += f"\n{body_text_text}"
                                else:
                                    description_text += body_text_text
                    if description_text != "":
                        output["description"] = description_text

            if bool(citation_dict):
                output["citation"] = citation_dict

            if "name" not in output:
                logger.info(f'no name for {output["identifier"]}')

            if "description" not in output:
                logger.info(f'no description for {output["identifier"]}')

            yield output
        logger.info(f"Articles not containing supplementary materials: {no_supplementary_material_count}")

    def update_cache(self):
        """If cache is not expired get the new records to add to the cache since last_updated"""

        last_updated = self.retreive_last_updated()

        new_files = self.get_zipped_files(last_updated)

        logger.info("Starting to download files")

        count = 0
        for zipped_filename in new_files:
            retry_count = 0
            while True:
                try:
                    logger.info("Downloading %s", zipped_filename)
                    wget.download("https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/xml/" + zipped_filename)
                    break
                except ContentTooShortError as e:
                    retry_count += 1
                    if retry_count > 10:
                        raise e
                    else:
                        logger.error("Error downloading file. %s. Retrying...", e)
                        logger.info("Try count: %s of 10", retry_count)
                        sleep(retry_count * 2)
            logger.info("Downloaded %s", zipped_filename)

            unzipped_file = tarfile.open(zipped_filename)
            logger.info("Unzipped %s", zipped_filename)

            unzipped_filepaths = unzipped_file.getnames()
            logger.info("Retrieved list of unzipped files in %s", zipped_filename)

            unzipped_file.extractall()
            logger.info("Extracted %s", zipped_filename)

            unzipped_file.close()
            logger.info("Closed %s", zipped_filename)

            os.remove(zipped_filename)
            logger.info("Deleted %s", zipped_filename)

            logger.info("Starting to yield files to cache")
            for xml_filepath in unzipped_filepaths:
                count += 1
                if count % 1000 == 0:
                    logger.info("Yielded %s files to cache", count)

                record_file = open(xml_filepath)
                # logger.info('Opened %s', xml_filepath)

                try:
                    record = ElementTree.parse(record_file)

                    root = record.getroot()
                    metadata_dict = self.get_metadata(root)
                    record_xml_string = ElementTree.tostring(root, encoding="unicode")
                    doc = {"metadata": metadata_dict, "xml": record_xml_string}

                    yield (xml_filepath.split("/")[1].replace(".xml", ""), json.dumps(doc))
                    # logger.info('Yielded %s', xml_filepath)

                    record_file.close()
                    # logger.info('Closed %s', xml_filepath)
                except ElementTree.ParseError as e:
                    logger.error("Error parsing %s", xml_filepath)
                    logger.error(e)
                    record_file.close()
                    os.remove(xml_filepath)
                    continue

            directories = set([x.split("/")[0] for x in unzipped_filepaths])
            [shutil.rmtree(x, ignore_errors=True) for x in directories]
            logger.info("Deleted directories: %s", directories)

        self.insert_last_updated()
