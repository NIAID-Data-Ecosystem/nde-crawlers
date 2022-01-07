# NIAID Data Ecosystem - Metadata Crawlers

This repository contains files to produce images that crawl various dataset
repositories to harvest their metadata.

Along with these also contain a BioThings Hub, which is used to ingest crawled
data, a BioThings Studio WebApp which simplifies using the Hub, and a Docker
Compose file that orchestrates multiple containers.

## Protocol

The crawler containers and the Hub exchange data by using a shared Docker
volume. Each crawler stores the metadata for a given source at 
`/data/<name_of_source>_crawled/data.ndjson` inside its container. Along with
`data.ndjson` should be `release.txt` which stores date and time of the start
of the metadata harvest, in ISO-8601 format.

When updating, the two files MUST not be overwritten. It should instead write
to temporary files and then rename the files when done. There are two reasons
for this. One, it makes sure that updates are atomic, either it fully succeeds
or not. Two, the Hub Dumper uses hard links to duplicate the data for uploading
so that file operations are fast and do not occupy unnecessary disk space,
overwriting the files will cause the data to be overwritten in the hub dumped
data archive.

For more details, see `NDEFileSystemDumper` and the code for harvesting data
from ImmPort.

## Current Issues
Due to the data transformation being incomplete. While data can be dumped and
uploaded in the Hub, it cannot be indexed to Elasticsearch as it conflicts with
the given mapping.