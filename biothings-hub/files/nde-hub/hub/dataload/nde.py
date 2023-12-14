import datetime
import inspect
import os
import shutil
import time

import orjson
from biothings.hub.dataload.dumper import BaseDumper
from biothings.hub.dataload.storage import IgnoreDuplicatedStorage
from biothings.hub.dataload.uploader import BaseSourceUploader
from config import CRAWLER_OUTPUT_DATA_ROOT, DATA_ARCHIVE_ROOT
from utils.utils import nde_upload_wrapper

__all__ = [
    "NDEFileSystemDumper",
    "NDESourceUploader",
]


class NDEFileSystemDumper(BaseDumper):
    """
    NDE Base Dumper

    This will create hard links for files from a given source
    """

    SCHEDULE = "0 5,17 * * *"  # 5am, 5pm

    @property
    def SRC_ROOT_FOLDER(self):
        if not self.SRC_NAME:
            raise RuntimeError
        return os.path.join(DATA_ARCHIVE_ROOT, self.SRC_NAME)

    def create_todump_list(self, force=False, **kwargs):
        """
        This gets called by method `dump`, to populate self.to_dump
        """
        # source location in same shared volume
        src_name = self.SRC_NAME
        release_filename = os.path.join(CRAWLER_OUTPUT_DATA_ROOT, f"{src_name}_crawled", "release.txt")
        if not os.path.exists(release_filename):
            # probably not crawled yet
            self.to_dump = []
            return
        data_filename = os.path.join(os.path.dirname(release_filename), "data.ndjson")
        # cat file and delay 5 seconds to sync release.txt issue with mnt from su07
        os.system("cat " + release_filename)
        time.sleep(5)

        # read crawled release string
        with open(release_filename, "r") as fh:
            next_release = fh.readline(63)
        # if a current dump exist, check release newer
        if not force and self.current_data_folder:
            existing_release_fn = os.path.join(self.current_data_folder, "release.txt")
            with open(existing_release_fn, "r") as fh:
                existing_release = fh.readline(63)
            ex_rl_dt = datetime.datetime.fromisoformat(existing_release.replace("Z", "+00:00"))
            nxt_rl_dt = datetime.datetime.fromisoformat(next_release.replace("Z", "+00:00"))
            # same release
            if ex_rl_dt == nxt_rl_dt:
                self.to_dump = []
                return
            elif ex_rl_dt > nxt_rl_dt:
                raise RuntimeError("Crawled release is older than dumped release")
            elif ex_rl_dt < nxt_rl_dt:
                pass
            else:
                raise RuntimeError("???")
        # either current dump does not exist or release is older
        self.release = next_release
        local_release_fn = os.path.join(self.new_data_folder, "release.txt")
        local_data_fn = os.path.join(self.new_data_folder, "data.ndjson")
        self.to_dump = [
            {"remote": release_filename, "local": os.path.abspath(local_release_fn)},
            {"remote": data_filename, "local": os.path.abspath(local_data_fn)},
        ]

    def remote_is_better(self, remotefile, localfile):
        """
        If there is a simple method to check whether remote is better
        """
        raise RuntimeError("remote_is_better is not available in NDEFileSystemDumper")

    def download(self, remotefile, localfile):
        """
        Create hard link
        """
        target_dir = os.path.dirname(localfile)
        os.makedirs(target_dir, exist_ok=True)
        try:
            os.link(remotefile, localfile)
        except OSError:
            # if link does not work, let's just download the file, e.g. when source/target is a mounted network folder
            shutil.copy(remotefile, localfile)

    @property
    def client(self):
        # overrides the parent class
        return None

    def prepare_client(self):
        raise RuntimeError("prepare_client method of NDEFileSystemDumper and its " "descendents must not be called")

    def release_client(self):
        # dump will always call this method so we have to allow it
        if inspect.stack()[1].function == "dump":
            return
        raise RuntimeError("release_client method of NDEFileSystemDumper and its " "descendents must not be called")

    def need_prepare(self):
        raise RuntimeError("need_prepare method of NDEFileSystemDumper and its " "descendents must not be called")


class NDESourceUploader(BaseSourceUploader):
    storage_class = IgnoreDuplicatedStorage

    @nde_upload_wrapper
    def load_data(self, data_folder):
        with open(os.path.join(data_folder, "data.ndjson"), "rb") as f:
            for line in f:
                doc = orjson.loads(line)
                yield doc

    @classmethod
    def get_mapping(cls):
        mapping = {
            "_meta": {
                "properties": {
                    "completeness": {
                        "properties": {
                            "augmented_recommended_ratio": {"type": "float"},
                            "augmented_required_ratio": {"type": "float"},
                            "recommended_max_score": {"type": "integer"},
                            "recommended_score": {"type": "integer"},
                            "recommended_score_ratio": {"type": "float"},
                            "required_max_score": {"type": "integer"},
                            "required_ratio": {"type": "float"},
                            "required_score": {"type": "integer"},
                            "total_max_score": {"type": "integer"},
                            "total_recommended_augmented": {"type": "integer"},
                            "total_required_augmented": {"type": "integer"},
                            "total_score": {"type": "integer"},
                            "weighted_score": {"type": "float"}
                        }
                    },
                    "recommended_augmented_fields": {"type": "keyword"},
                    "required_augmented_fields": {"type": "keyword"}
                }
            },
            "@context": {
                "type": "text",
                "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
            },
            "@id": {"type": "keyword"},
            "@type": {"type": "keyword", "copy_to": ["all"]},
            "abstract": {"type": "text", "analyzer": "nde_analyzer", "copy_to": ["all"]},
            "all": {
                "type": "text",
                "analyzer": "nde_analyzer",
                "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
            },
            "aggregateRating": {
                "properties": {
                    "@type": {"type": "text"},
                    "ratingCount": {"type": "unsigned_long"},
                    "ratingValue": {"type": "double"},
                    "reviewAspect": {"type": "text"},
                }
            },
            "alternateName": {"type": "text", "copy_to": ["all"]},
            "applicationCategory": {
                "type": "keyword",
                "normalizer": "keyword_lowercase_normalizer",
                "copy_to": ["all"],
            },
            "applicationSubCategory": {
                "properties": {
                    "@type": {"type": "text"},
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "keyword", "copy_to": ["all"]},
                    "url": {"type": "keyword"},
                }
            },
            "applicationSuite": {"type": "text", "analyzer": "nde_analyzer", "copy_to": ["all"]},
            "author": {
                "properties": {
                    "@type": {"type": "text"},
                    "affiliation": {
                        "properties": {
                            "@type": {"type": "keyword", "copy_to": ["all"]},
                            "name": {"type": "keyword", "copy_to": ["all"]},
                            "sameAs": {"type": "keyword", "copy_to": ["all"]},
                        }
                    },
                    "email": {"type": "text", "copy_to": ["all"]},
                    "familyName": {"type": "text", "copy_to": ["all"]},
                    "givenName": {"type": "text", "copy_to": ["all"]},
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "text", "copy_to": ["all"]},
                    "role": {"type": "keyword"},
                    "title": {"type": "text"},
                    "url": {"type": "keyword"},
                }
            },
            "availableOnDevice": {"type": "text"},
            "citation": {
                "properties": {
                    "@type": {"type": "keyword"},
                    "abstract": {"type": "text"},
                    "author": {
                        "properties": {
                            "@type": {"type": "text"},
                            "familyName": {"type": "text", "copy_to": ["all"]},
                            "givenName": {"type": "text", "copy_to": ["all"]},
                            "name": {"type": "text", "copy_to": ["all"]},
                        }
                    },
                    "citation": {"type": "text"},
                    "datePublished": {"type": "date"},
                    "description": {"type": "text"},
                    "doi": {"type": "keyword", "copy_to": ["all"]},
                    "identifier": {"type": "keyword", "copy_to": ["all"]},
                    "issueNumber": {"type": "text"},
                    "journalName": {"type": "keyword", "copy_to": ["all"]},
                    "journalNameAbbrev": {"type": "keyword", "copy_to": ["all"]},
                    "name": {"type": "text", "analyzer": "nde_analyzer", "copy_to": ["all"]},
                    "pagination": {"type": "text"},
                    "pmid": {"type": "text", "copy_to": ["all"]},
                    "url": {"type": "keyword"},
                    "volumeNumber": {"type": "text"},
                }
            },
            "citedBy": {
                "properties": {
                    "@type": {"type": "keyword"},
                    "abstract": {"type": "text"},
                    "author": {
                        "properties": {
                            "@type": {"type": "text"},
                            "familyName": {"type": "text", "copy_to": ["all"]},
                            "givenName": {"type": "text", "copy_to": ["all"]},
                            "name": {"type": "text", "copy_to": ["all"]},
                        }
                    },
                    "citation": {"type": "text"},
                    "datePublished": {"type": "date"},
                    "description": {"type": "text"},
                    "doi": {"type": "keyword", "copy_to": ["all"]},
                    "identifier": {"type": "keyword", "copy_to": ["all"]},
                    "issueNumber": {"type": "text"},
                    "journalName": {"type": "keyword", "copy_to": ["all"]},
                    "journalNameAbbrev": {"type": "keyword", "copy_to": ["all"]},
                    "name": {"type": "text", "analyzer": "nde_analyzer", "copy_to": ["all"]},
                    "pagination": {"type": "text"},
                    "pmid": {"type": "text", "copy_to": ["all"]},
                    "url": {"type": "keyword"},
                    "volumeNumber": {"type": "text"},
                }
            },
            "codeRepository": {
                "type": "text",
                "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
            },
            "conditionsOfAccess": {"type": "keyword"},
            "contentSize": {"type": "text"},
            "contentUrl": {
                "type": "text",
                "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
            },
            "contributor": {
                "properties": {
                    "@id": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "@type": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "affiliation": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "email": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "name": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                }
            },
            "creator": {
                "properties": {
                    "@type": {"type": "text"},
                    "affiliation": {
                        "properties": {
                            "@type": {"type": "keyword", "copy_to": ["all"]},
                            "name": {"type": "keyword", "copy_to": ["all"]},
                            "sameAs": {"type": "keyword", "copy_to": ["all"]},
                        }
                    },
                    "email": {"type": "text", "copy_to": ["all"]},
                    "familyName": {"type": "text", "copy_to": ["all"]},
                    "givenName": {"type": "text", "copy_to": ["all"]},
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "text", "copy_to": ["all"]},
                    "role": {"type": "keyword"},
                    "title": {"type": "text"},
                    "url": {"type": "keyword"},
                }
            },
            "curatedBy": {
                "properties": {
                    "@type": {"type": "text"},
                    "name": {"type": "keyword", "copy_to": ["all"]},
                    "url": {"type": "text"},
                    "versionDate": {"type": "date"},
                }
            },
            "date": {"type": "date"},
            "dateCreated": {"type": "date"},
            "dateModified": {"type": "date"},
            "datePublished": {"type": "date"},
            "description": {"type": "text", "analyzer": "nde_analyzer", "copy_to": ["all"]},
            "discussionUrl": {"type": "text", "copy_to": ["all"]},
            "distribution": {
                "properties": {
                    "@id": {"type": "keyword"},
                    "@type": {"type": "keyword"},
                    "contentUrl": {"type": "text"},
                    "dateCreated": {"type": "date"},
                    "dateModified": {"type": "date"},
                    "datePublished": {"type": "date"},
                    "description": {"type": "text", "analyzer": "nde_analyzer"},
                    "encodingFormat": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "contentSize": {"type": "text"},
                    "name": {"type": "keyword"},
                }
            },
            "doi": {"type": "text", "copy_to": ["all"], "fields": {"keyword": {"type": "keyword"}}},
            "downloadUrl": {
                "properties": {
                    "name": {"type": "text", "copy_to": ["all"]},
                }
            },
            "featureList": {"type": "text", "copy_to": ["all"]},
            "funding": {
                "properties": {
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "url": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "keyword", "normalizer": "keyword_lowercase_normalizer", "copy_to": ["all"]},
                    "description": {"type": "text", "copy_to": ["all"]},
                    "endDate": {"type": "date"},
                    "startDate": {"type": "date"},
                    "keywords": {"type": "keyword", "normalizer": "keyword_lowercase_normalizer", "copy_to": ["all"]},
                    "isBasedOn": {
                        "properties": {
                            "identifier": {"type": "text", "copy_to": ["all"]},
                        }
                    },
                    "funder": {
                        "properties": {
                            "alternateName": {"type": "keyword", "copy_to": ["all"]},
                            "class": {"type": "keyword"},
                            "description": {"type": "text", "copy_to": ["all"]},
                            "identifier": {"type": "text", "copy_to": ["all"]},
                            "name": {
                                "type": "keyword",
                                "normalizer": "keyword_lowercase_normalizer",
                                "copy_to": ["all"],
                            },
                            "parentOrganization": {"type": "keyword", "copy_to": ["all"]},
                            "role": {"type": "keyword"},
                            "url": {"type": "text", "copy_to": ["all"]},
                            "employee": {
                                "properties": {
                                    "givenName": {"type": "text", "copy_to": ["all"]},
                                    "familyName": {"type": "text", "copy_to": ["all"]},
                                    "name": {
                                        "type": "keyword",
                                        "normalizer": "keyword_lowercase_normalizer",
                                        "copy_to": ["all"],
                                    },
                                }
                            },
                        },
                    },
                }
            },
            "hasPart": {
                "properties": {
                    "@id": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "@type": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "additionalType": {
                        "properties": {
                            "name": {"type": "keyword", "copy_to": ["all"]},
                            "url": {"type": "text"},
                        }
                    },
                    "encodingFormat": {"type": "text"},
                    "name": {"type": "text", "copy_to": ["all"]},
                    "url": {"type": "text"},
                }
            },
            "healthCondition": {
                "properties": {
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "inDefinedTermSet": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "keyword", "normalizer": "keyword_lowercase_normalizer", "copy_to": ["all"]},
                    "url": {"type": "text", "copy_to": ["all"]},
                    "alternateName": {
                        "type": "keyword",
                        "normalizer": "keyword_lowercase_normalizer",
                        "copy_to": ["all"],
                    },
                    "originalName": {
                        "type": "keyword",
                        "normalizer": "keyword_lowercase_normalizer",
                        "copy_to": ["all"],
                    },
                    "isCurated": {"type": "boolean"},
                    "curatedBy": {
                        "properties": {
                            "name": {"type": "keyword", "copy_to": ["all"]},
                            "url": {"type": "text", "copy_to": ["all"]},
                            "dateModified": {"type": "date"},
                        }
                    },
                }
            },
            "identifier": {"type": "text", "copy_to": ["all"]},
            "input": {
                "properties": {
                    "@type": {"type": "keyword"},
                    "description": {"type": "text"},
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "text", "copy_to": ["all"]},
                    "encodingFormat": {"type": "text", "copy_to": ["all"]},
                },
            },
            "includedInDataCatalog": {
                "properties": {
                    "@type": {"type": "text"},
                    "name": {"type": "keyword", "copy_to": ["all"]},
                    "url": {"type": "text"},
                    "versionDate": {"type": "date"},
                }
            },
            "inLanguage": {
                "properties": {
                    "@type": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "alternateName": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "name": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                }
            },
            "infectiousAgent": {
                "properties": {
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "inDefinedTermSet": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "keyword", "normalizer": "keyword_lowercase_normalizer", "copy_to": ["all"]},
                    "url": {"type": "text", "copy_to": ["all"]},
                    "alternateName": {
                        "type": "keyword",
                        "normalizer": "keyword_lowercase_normalizer",
                        "copy_to": ["all"],
                    },
                    "originalName": {
                        "type": "keyword",
                        "normalizer": "keyword_lowercase_normalizer",
                        "copy_to": ["all"],
                    },
                    "commonName": {"type": "keyword", "normalizer": "keyword_lowercase_normalizer", "copy_to": ["all"]},
                    "displayName": {
                        "type": "keyword",
                        "normalizer": "keyword_lowercase_normalizer",
                        "copy_to": ["all"],
                    },
                    "isCurated": {"type": "boolean"},
                    "curatedBy": {
                        "properties": {
                            "name": {"type": "keyword", "copy_to": ["all"]},
                            "url": {"type": "text", "copy_to": ["all"]},
                            "dateModified": {"type": "date"},
                        }
                    },
                }
            },
            "interactionStatistic": {
                "properties": {
                    "@type": {"type": "text"},
                    "interactionType": {"type": "text"},
                    "userInteractionCount": {"type": "unsigned_long"},
                }
            },
            "isAccessibleForFree": {"type": "boolean"},
            "isBasedOn": {
                "properties": {
                    "@type": {"type": "keyword"},
                    "abstract": {"type": "text"},
                    "additionalType": {
                        "properties": {
                            "name": {"type": "keyword", "copy_to": ["all"]},
                            "url": {"type": "text"},
                        }
                    },
                    "citation": {"type": "text"},
                    "datePublished": {"type": "date"},
                    "description": {"type": "text", "analyzer": "nde_analyzer"},
                    "doi": {"type": "text", "copy_to": ["all"]},
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "text", "analyzer": "nde_analyzer"},
                    "pmid": {"type": "text", "copy_to": ["all"]},
                    "url": {"type": "text"},
                }
            },
            "isBasisFor": {
                "properties": {
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "text", "copy_to": ["all"]},
                    "url": {"type": "text"},
                }
            },
            "isPartOf": {
                "properties": {
                    "@id": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "@type": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "alternateName": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "text", "copy_to": ["all"]},
                    "identifier": {"type": "keyword", "copy_to": ["all"]},
                    "url": {"type": "text"},
                }
            },
            "isRelatedTo": {
                "properties": {
                    "@type": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "citation": {
                        "properties": {
                            "pmid": {"type": "text", "copy_to": ["all"]},
                            "url": {"type": "text"},
                        }
                    },
                    "name": {"type": "text", "analyzer": "nde_analyzer", "copy_to": ["all"]},
                    "identifier": {"type": "keyword", "copy_to": ["all"]},
                    "includedInDataCatalog": {
                        "properties": {
                            "@type": {"type": "text"},
                            "name": {"type": "keyword", "copy_to": ["all"]},
                            "url": {"type": "text"},
                            "versionDate": {"type": "date"},
                        }
                    },
                    "hasPart": {
                        "properties": {
                            "@type": {"type": "text"},
                            "identifier": {"type": "keyword", "copy_to": ["all"]},
                        }
                    },
                    "relationship": {"type": "text", "copy_to": ["all"]},
                    "url": {"type": "text"},
                }
            },
            "isSimilarTo": {
                "properties": {
                    "@type": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "name": {"type": "text", "analyzer": "nde_analyzer", "copy_to": ["all"]},
                    "identifier": {"type": "keyword", "copy_to": ["all"]},
                    "includedInDataCatalog": {
                        "properties": {
                            "@type": {"type": "text"},
                            "name": {"type": "keyword", "copy_to": ["all"]},
                            "url": {"type": "text"},
                            "versionDate": {"type": "date"},
                        }
                    },
                    "relationship": {"type": "keyword", "copy_to": ["all"]},
                }
            },
            "keywords": {"type": "keyword", "normalizer": "keyword_lowercase_normalizer", "copy_to": ["all"]},
            "license": {"type": "text"},
            "measurementTechnique": {
                "properties": {
                    "description": {"type": "text"},
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "keyword", "normalizer": "keyword_lowercase_normalizer", "copy_to": ["all"]},
                    "url": {"type": "text", "copy_to": ["all"]},
                }
            },
            "mainEntityOfPage": {"type": "text"},
            "metadata_score": {"type": "float"},
            "name": {
                "type": "text",
                "analyzer": "nde_analyzer",
                "copy_to": ["all"],
                "fields": {
                    "phrase_suggester": {"type": "text", "analyzer": "phrase_suggester"},
                    "raw": {"type": "keyword"},
                },
            },
            "nctid": {"type": "keyword", "copy_to": ["all"]},
            "output": {
                "properties": {
                    "@type": {"type": "keyword"},
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "text", "analyzer": "nde_analyzer", "copy_to": ["all"]},
                    "encodingFormat": {"type": "text", "copy_to": ["all"]},
                },
            },
            "processorRequirements": {"type": "text", "copy_to": ["all"]},
            "programmingLanguage": {
                "type": "keyword",
                "normalizer": "keyword_lowercase_normalizer",
                "copy_to": ["all"],
            },
            "relationship": {"type": "text", "copy_to": ["all"]},
            "sameAs": {
                "type": "text",
                "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
            },
            "sdPublisher": {
                "properties": {
                    "@type": {"type": "keyword", "copy_to": ["all"]},
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "text", "copy_to": ["all"]},
                    "url": {"type": "keyword"},
                }
            },
            "softwareAddOn": {
                "properties": {
                    "identifier": {"type": "text", "copy_to": ["all"]},
                }
            },
            "softwareHelp": {
                "properties": {
                    "name": {"type": "text", "copy_to": ["all"]},
                    "url": {"type": "text", "copy_to": ["all"]},
                }
            },
            "softwareRequirements": {"type": "text", "analyzer": "nde_analyzer", "copy_to": ["all"]},
            "softwareVersion": {"type": "text", "copy_to": ["all"]},
            "sourceInfo": {
                "properties": {
                    "name": {"type": "text"},
                    "description": {"type": "text"},
                    "schema": {"properties": {"field": {"type": "text"}}},
                    "url": {"type": "text"},
                    "identifier": {"type": "text"},
                }
            },
            "spatialCoverage": {
                "properties": {
                    "@type": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "geo": {
                        "properties": {
                            "@type": {
                                "type": "text",
                                "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                            },
                            "latitude": {"type": "float"},
                            "longitude": {"type": "float"},
                        }
                    },
                    "name": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "identifier": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                }
            },
            "species": {
                "properties": {
                    "additionalType": {
                        "properties": {
                            "name": {"type": "keyword", "copy_to": ["all"]},
                            "url": {"type": "text"},
                        }
                    },
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "inDefinedTermSet": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "keyword", "normalizer": "keyword_lowercase_normalizer", "copy_to": ["all"]},
                    "url": {"type": "text", "copy_to": ["all"]},
                    "alternateName": {
                        "type": "keyword",
                        "normalizer": "keyword_lowercase_normalizer",
                        "copy_to": ["all"],
                    },
                    "originalName": {
                        "type": "keyword",
                        "normalizer": "keyword_lowercase_normalizer",
                        "copy_to": ["all"],
                    },
                    "commonName": {"type": "keyword", "normalizer": "keyword_lowercase_normalizer", "copy_to": ["all"]},
                    "displayName": {
                        "type": "keyword",
                        "normalizer": "keyword_lowercase_normalizer",
                        "copy_to": ["all"],
                    },
                    "isCurated": {"type": "boolean"},
                    "curatedBy": {
                        "properties": {
                            "name": {"type": "keyword", "copy_to": ["all"]},
                            "url": {"type": "text", "copy_to": ["all"]},
                            "dateModified": {"type": "date"},
                        }
                    },
                }
            },
            "temporalCoverage": {
                "properties": {
                    "temporalInterval": {
                        "properties": {
                            "@type": {"type": "text"},
                            "duration": {"type": "text", "copy_to": ["all"]},
                            "endDate": {"type": "date"},
                            "name": {"type": "text"},
                            "startDate": {"type": "date"},
                        }
                    }
                }
            },
            "thumbnailUrl": {
                "type": "text",
                "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
            },
            "topicCategory": {
                "properties": {
                    "description": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "keyword", "copy_to": ["all"]},
                    "url": {"type": "text", "copy_to": ["all"]},
                    "curatedBy": {
                        "properties": {
                            "name": {"type": "text", "copy_to": ["all"]},
                            "url": {"type": "text", "copy_to": ["all"]},
                        }
                    },
                }
            },
            "url": {"type": "text", "copy_to": ["all"]},
            "usageInfo": {
                "properties": {
                    "description": {"type": "text", "analyzer": "nde_analyzer", "copy_to": ["all"]},
                    "name": {"type": "text", "analyzer": "nde_analyzer", "copy_to": ["all"]},
                    "url": {"type": "text", "copy_to": ["all"]},
                }
            },
            "variableMeasured": {"type": "keyword", "normalizer": "keyword_lowercase_normalizer", "copy_to": ["all"]},
            "version": {
                "type": "text",
                "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
            },
        }

        return mapping
