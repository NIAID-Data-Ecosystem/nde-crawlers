import datetime
import inspect
import os
import shutil

import orjson
from biothings.hub.dataload.dumper import BaseDumper
from biothings.hub.dataload.storage import IgnoreDuplicatedStorage
from biothings.hub.dataload.uploader import BaseSourceUploader

from utils.date import add_date
from config import DATA_ARCHIVE_ROOT, CRAWLER_OUTPUT_DATA_ROOT


__all__ = [
    'NDEFileSystemDumper',
    'NDESourceUploader',
]


class NDEFileSystemDumper(BaseDumper):
    """
    NDE Base Dumper

    This will create hard links for files from a given source
    """
    SCHEDULE = "0 * * * *"  # hourly

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
        release_filename = os.path.join(
            CRAWLER_OUTPUT_DATA_ROOT, f'{src_name}_crawled', 'release.txt'
        )
        if not os.path.exists(release_filename):
            # probably not crawled yet
            self.to_dump = []
            return
        data_filename = os.path.join(
            os.path.dirname(release_filename), 'data.ndjson'
        )
        # read crawled release string
        with open(release_filename, 'r') as fh:
            next_release = fh.readline(63)
        # if a current dump exist, check release newer
        if not force and self.current_data_folder:
            existing_release_fn = os.path.join(
                self.current_data_folder, 'release.txt'
            )
            with open(existing_release_fn, 'r') as fh:
                existing_release = fh.readline(63)
            ex_rl_dt = datetime.datetime.fromisoformat(
                existing_release.replace('Z', '+00:00')
            )
            nxt_rl_dt = datetime.datetime.fromisoformat(
                next_release.replace('Z', '+00:00')
            )
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
        local_release_fn = os.path.join(
            self.new_data_folder, 'release.txt'
        )
        local_data_fn = os.path.join(
            self.new_data_folder, 'data.ndjson'
        )
        self.to_dump = [
            {
                'remote': release_filename,
                'local': os.path.abspath(local_release_fn)
            },
            {
                'remote': data_filename,
                'local': os.path.abspath(local_data_fn)
            }
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
        raise RuntimeError("prepare_client method of NDEFileSystemDumper and its "
                           "descendents must not be called")

    def release_client(self):
        # dump will always call this method so we have to allow it
        if inspect.stack()[1].function == 'dump':
            return
        raise RuntimeError("release_client method of NDEFileSystemDumper and its "
                           "descendents must not be called")

    def need_prepare(self):
        raise RuntimeError("need_prepare method of NDEFileSystemDumper and its "
                           "descendents must not be called")


class NDESourceUploader(BaseSourceUploader):
    storage_class = IgnoreDuplicatedStorage

    def load_data(self, data_folder):
        with open(os.path.join(data_folder, 'data.ndjson'), 'rb') as f:
            for line in f:
                doc = orjson.loads(line)
                #add date transformation here and have the most recent date of all of the dates
                doc = add_date(doc)
                yield doc

    @classmethod
    def get_mapping(cls):
        mapping = {
            "@context": {
                "type": "text",
                "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
            },
            "@id": {"type": "keyword"},
            "@type": {"type": "keyword", "copy_to": ["all"]},
            "about": {
                "properties": {
                    "@id": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "@type": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                }
            },
            "abstract": {"type": "text", "copy_to": ["all"]},
            "all": {
                "type": "text",
                "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
            },
            "alternateName": {"type": "text", "copy_to": ["all"]},
            "armGroup": {
                "properties": {
                    "@type": {"type": "keyword"},
                    "description": {"type": "text", "copy_to": ["all"]},
                    "intervention": {
                        "properties": {
                            "@type": {"type": "keyword"},
                            "category": {"type": "keyword", "copy_to": ["all"]},
                            "description": {"type": "text", "copy_to": ["all"]},
                            "name": {"type": "keyword", "copy_to": ["all"]},
                        }
                    },
                    "name": {"type": "keyword", "copy_to": ["all"]},
                    "role": {"type": "keyword", "copy_to": ["all"]},
                }
            },
            "author": {
                "properties": {
                    "@type": {"type": "text"},
                    "affiliation": {
                        "properties": {"name": {"type": "keyword", "copy_to": ["all"]}}
                    },
                    "familyName": {"type": "text", "copy_to": ["all"]},
                    "givenName": {"type": "text", "copy_to": ["all"]},
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "text", "copy_to": ["all"]},
                    "role": {"type": "keyword"},
                    "title": {"type": "text"},
                    "url": {"type": "keyword"},
                }
            },
            "citation": {
                "properties": {
                    "@type": {"type": "keyword"},
                    "author": {
                        "properties": {
                            "@type": {"type": "text"},
                            "familyName": {"type": "text", "copy_to": ["all"]},
                            "givenName" : {"type": "text", "copy_to": ["all"]},
                            "name": {"type": "text", "copy_to": ["all"]},
                        }
                    },
                    "datePublished": {"type": "date"},
                    "doi": {"type": "keyword", "copy_to": ["all"]},
                    "identifier": {"type": "keyword", "copy_to": ["all"]},
                    "issueNumber": {"type": "text"},
                    "journalName": {"type": "keyword", "copy_to": ["all"]},
                    "journalNameAbbrev": {"type": "keyword", "copy_to": ["all"]},
                    "name": {"type": "text", "copy_to": ["all"]},
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
                    "citation": {"type": "text"},
                    "datePublished": {"type": "text"},
                    "description": {"type": "text"},
                    "doi": {"type": "text", "copy_to": ["all"]},
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "text"},
                    "pmid": {"type": "text", "copy_to": ["all"]},
                    "url": {"type": "text"},
                }
            },
            "codeRepository": {
                "type": "text",
                "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
            },
            "conditionsOfAccess": {"type": "text"},
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
                    "name": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                }
            },
            "correction": {"type": "text"},
            "creator": {
                "properties": {
                    "@id": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "@type": {"type": "text"},
                    "affiliation": {
                        "properties": {"name": {"type": "keyword", "copy_to": ["all"]}}
                    },
                    "familyName": {"type": "text", "copy_to": ["all"]},
                    "givenName": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "text", "copy_to": ["all"]},
                    "role": {"type": "keyword"},
                    "title": {"type": "text"},
                }
            },
            # TODO: AFTER FIXING ZENODO DELETE THIS MAPPING REPLACED BY includedInDataCatalog
            "curatedBy": {
                "properties": {
                    "@type": {"type": "text"},
                    "identifier": {"type": "keyword"},
                    "name": {"type": "keyword", "copy_to": ["all"]},
                    "url": {"type": "text"},
                    "versionDate": {"type": "keyword"},
                }
            },
            "date": {"type": "date"},
            "dateCreated": {"type": "keyword"},
            "dateModified": {"type": "keyword"},
            "datePublished": {"type": "keyword"},
            "description": {"type": "text", "copy_to": ["all"]},
            "distribution": {
                "properties": {
                    "@id": {"type": "keyword"},
                    "@type": {"type": "keyword"},
                    "contentUrl": {"type": "text"},
                    "dateCreated": {"type": "keyword"},
                    "dateModified": {"type": "keyword"},
                    "datePublished": {"type": "keyword"},
                    "description": {"type": "text"},
                    "encodingFormat": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "name": {"type": "keyword"},
                }
            },
            "doi": {"type": "text", "copy_to": ["all"]},
            "duration": {"type": "text"},
            "eligibilityCriteria": {
                "properties": {
                    "@type": {"type": "keyword"},
                    "exclusionCriteria": {"type": "text", "copy_to": ["all"]},
                    "gender": {"type": "text"},
                    "healthyVolunteers": {"type": "boolean"},
                    "inclusionCriteria": {"type": "text", "copy_to": ["all"]},
                    "maximumAge": {"type": "text"},
                    "minimumAge": {"type": "text"},
                    "stdAge": {"type": "text", "copy_to": ["all"]},
                }
            },
            "funding": {
                "properties": {
                    "description": {"type": "text", "copy_to": ["all"]},
                    "funder": {
                        "properties": {
                            "alternateName": {"type": "keyword", "copy_to": ["all"]},
                            "class": {"type": "keyword"},
                            "name": {"type": "keyword", "copy_to": ["all"]},
                            "parentOrganization": {"type": "keyword", "copy_to": ["all"]},
                            "role": {"type": "keyword"},
                            "url": {"type": "text", "copy_to": ["all"]},
                        }
                    },
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "url": {"type": "text", "copy_to": ["all"]},
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
                }
            },
            "hasResults": {"type": "boolean"},
            "headline": {
                "type": "text",
                "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
            },
            "healthCondition": {"type": "text", "copy_to": ["all"]},
            "identifier": {"type": "text", "copy_to": ["all"]},
            "identifierSource": {"type": "keyword", "copy_to": ["all"]},
            "image": {
                "type": "text",
                "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
            },
            "includedInDataCatalog": {
                "properties": {
                    "@type": {"type": "text"},
                    "name": {"type": "keyword", "copy_to": ["all"]},
                    "url": {"type": "text"},
                    "versionDate": {"type": "keyword"},
                }
            },
            "inComplianceWith": {"type": "text"},
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
                     "name": {"type": "keyword", "copy_to": ["all"]},
                     "url": {"type": "text", "copy_to": ["all"]},
                }
            },
            "infectiousDisease": {
                "properties": {
                     "name": {"type": "keyword", "copy_to": ["all"]},
                     "url": {"type": "text", "copy_to": ["all"]},
                }
            },
            "instrument": {"type": "text", "copy_to": ["all"]},
            "interventionText": {"type": "text", "copy_to": ["all"]},
            "interventions": {
                "properties": {
                    "@type": {"type": "keyword"},
                    "category": {"type": "keyword", "copy_to": ["all"]},
                    "description": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "keyword", "copy_to": ["all"]},
                }
            },
            "isBasedOn": {
                "properties": {
                    "@type": {"type": "keyword"},
                    "abstract": {"type": "text"},
                    "citation": {"type": "text"},
                    "datePublished": {"type": "text"},
                    "description": {"type": "text"},
                    "doi": {"type": "text", "copy_to": ["all"]},
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "text"},
                    "pmid": {"type": "text", "copy_to": ["all"]},
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
                }
            },
            "issueNumber": {"type": "text"},
            "journalName": {"type": "keyword", "copy_to": ["all"]},
            "journalNameAbbrev": {"type": "keyword", "copy_to": ["all"]},
            "keywords": {"type": "keyword", "copy_to": ["all"]},
            "license": {"type": "text"},
            "material": {"type": "text", "copy_to": ["all"]},
            "measurementParameter": {"properties": {"resolution": {"type": "keyword"}}},
            "measurementTechnique": {                
                "properties": {
                     "name": {"type": "keyword", "copy_to": ["all"]},
                     "url": {"type": "text", "copy_to": ["all"]},
                }
            },
            "name": {"type": "keyword", "copy_to": ["all"]},
            "outcome": {
                "properties": {
                    "@type": {"type": "keyword"},
                    "outcomeMeasure": {"type": "text", "copy_to": ["all"]},
                    "outcomeTimeFrame": {"type": "text", "copy_to": ["all"]},
                    "outcomeType": {"type": "keyword", "copy_to": ["all"]},
                }
            },
            "pmid": {"type": "integer", "copy_to": ["all"]},
            "protocolCategory": {"type": "keyword", "copy_to": ["all"]},
            "protocolSetting": {"type": "keyword", "copy_to": ["all"]},
            "publicationType": {"type": "keyword", "copy_to": ["all"]},
            "relatedTo": {
                "properties": {
                    "@type": {"type": "keyword"},
                    "abstract": {"type": "text"},
                    "citation": {"type": "text"},
                    "datePublished": {"type": "text"},
                    "description": {"type": "text"},
                    "doi": {"type": "text", "copy_to": ["all"]},
                    "identifier": {"type": "text", "copy_to": ["all"]},
                    "name": {"type": "text", "copy_to": ["all"]},
                    "pmid": {"type": "text"},
                    "url": {"type": "text"},
                }
            },
            "sameAs": {
                "type": "text",
                "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
            },
            "spatial": {
                "properties": {
                    "@type": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "geo": {
                        "properties": {
                            "@type": {
                                "type": "text",
                                "fields": {
                                    "keyword": {"type": "keyword", "ignore_above": 256}
                                },
                            },
                            "latitude": {"type": "float"},
                            "longitude": {"type": "float"},
                        }
                    },
                    "name": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                }
            },
            "species": {
                "properties": {
                     "name": {"type": "keyword", "copy_to": ["all"]},
                     "url": {"type": "text", "copy_to": ["all"]},
                }
            },
            "sponsor": {
                "properties": {
                    "@type": {"type": "keyword"},
                    "class": {"type": "keyword", "copy_to": ["all"]},
                    "name": {"type": "keyword", "copy_to": ["all"]},
                    "role": {"type": "keyword"},
                }
            },
            "studyDesign": {
                "properties": {
                    "@type": {"type": "keyword"},
                    "designAllocation": {"type": "keyword", "copy_to": ["all"]},
                    "designModel": {"type": "keyword", "copy_to": ["all"]},
                    "designPrimaryPurpose": {"type": "keyword", "copy_to": ["all"]},
                    "designWhoMasked": {"type": "text", "copy_to": ["all"]},
                    "phase": {"type": "keyword", "copy_to": ["all"]},
                    "phaseNumber": {"type": "half_float"},
                    "studyDesignText": {"type": "text", "copy_to": ["all"]},
                    "studyType": {"type": "keyword", "copy_to": ["all"]},
                }
            },
            "studyEvent": {
                "properties": {
                    "@type": {"type": "keyword"},
                    "studyEventDate": {"type": "text"},
                    "studyEventDateType": {"type": "text"},
                    "studyEventType": {"type": "text"},
                }
            },
            "studyLocation": {
                "properties": {
                    "@type": {"type": "keyword"},
                    "name": {"type": "keyword", "copy_to": ["all"]},
                    "studyLocationCity": {"type": "keyword", "copy_to": ["all"]},
                    "studyLocationCountry": {"type": "keyword", "copy_to": ["all"]},
                    "studyLocationState": {"type": "keyword", "copy_to": ["all"]},
                    "studyLocationStatus": {"type": "keyword", "copy_to": ["all"]},
                }
            },
            "studyStatus": {
                "properties": {
                    "@type": {"type": "keyword"},
                    "enrollmentCount": {"type": "integer"},
                    "enrollmentType": {"type": "text"},
                    "status": {"type": "keyword", "copy_to": ["all"]},
                    "statusDate": {"type": "text"},
                    "statusExpanded": {"type": "boolean"},
                    "whyStopped": {"type": "text", "copy_to": ["all"]},
                }
            },
            "thumbnailUrl": {
                "type": "text",
                "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
            },
            "topicCategory": {"type": "keyword", "copy_to": ["all"]},
            "url": {"type": "text", "copy_to": ["all"]},
            "usedToGenerate": {"type": "text"},
            "variableMeasured": {"type": "keyword", "copy_to": ["all"]},
            "version": {
                "type": "text",
                "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
            },
            "volumeNumber": {"type": "text"},
            "workFeatured": {
                "properties": {
                    "@type": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "alternateName": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "location": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "name": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                    "url": {
                        "type": "text",
                        "fields": {"keyword": {"type": "keyword", "ignore_above": 256}},
                    },
                }
            },
        }

        return mapping
