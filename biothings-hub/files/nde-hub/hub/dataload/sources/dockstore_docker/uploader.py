from hub.dataload.nde import NDESourceUploader
from utils.csv_helper import get_source_data

# Example __metadata__ dictionary:
# <SOURCE_NAME> = https://api.data.niaid.nih.gov/v1/metadata
# __metadata__ = {
# "src_meta": {
# 'description': 'A short description of what the source offers, usually found on the source's about page',
# 'name': 'The full source name, Ex. Mendeley Data (not mendeley)',
# 'identifier':'includedInDataCatalog.name value',
# 'schema': 'A dict where the key is the source's metadata variable and the value is our transformation. Ex: {"summary":"description"},
# 'url': 'The source's URL homepage',
# }
# }


class DockStoreDockerUploader(NDESourceUploader):
    name = "dockstore_docker"
    __metadata__ = {
        "src_meta": {
            "sourceInfo": {
                "description": "Dockstore provides a place where users can share tools encapsulated in Docker and described with the Common Workflow Language(CWL) or Workflow Description Language(WDL), workflow languages used by members of and APIs created by the GA4GH Cloud Work Stream. This enables scientists, for example, to share analytical tools in a way that makes them machine readable and runnable in a variety of environments. While the Dockstore is focused on serving researchers in the biosciences, the combination of Docker + CWL/WDL can be used by anyone to describe the tools and services in their Docker images in a standardized, machine-readable way. Dockstore also attempts to work with new and alternative languages/standards such as Nextflow as popular challengers to CWL and WDL. While contributing work on the GA4GH Tool Registry standard as a way of sharing data with workflow platforms and partners.",
                "identifier": "Dockstore",
                "name": "Dockstore",
                "schema": get_source_data(name),
                "url": "https://dockstore.org/",
            }
        }
    }
