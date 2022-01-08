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

### How does this all work? An example

Let's use ImmPort as an example. There is a lot going on in the `Dockerfile`,
but deep down, aside from copying the code, it sets up the environment specific
to ImmPort crawling, which consists of things like the Python dependencies,
Google Chrome, and Chromedriver. Then it sets up the BioThings user which is
common across all containers so it won't run into permission issues when
creating files and directories. Then it installs `cron`, so that we can have
tasks run periodically.

The container starts by running `cron` as its only foreground task. It reads
`/etc/crontab` and `crontab` file asks cron to run our shell script
`run-spider.sh` on startup and also periodically (as root, explained later).

Our script starts the actual `spider.py` by running the `scrapy` command which
was installed as a dependency. However it does this with three twists: 1) it 
obtains a lock file `/crawler.lock` so that only one instance would be running
at once. 2) it redirects the standard output and standard error to that of the
PID 1 process (`cron`) so that any output of our script propagates to the output
streams that Docker can see, this is done as root because cron has to be run
as root and therefore only root user can write to the output of `cron` 3) run
the actual scrapy spider as the `biothings` user.

Our spider performs the crawling and finally produces 
`/data/immport_crawled/data.ndjson` and `/data/immport_crawled/release.txt`.
Here `/data` is a shared volume that all containers should have access to.
After a while, as the all descendants of `NDEFileSystemDumper` runs every hour,
the related dumper in
`biothings-hub/files/nde-hub/hub/dataload/sources/immport/dumper.py` will pick
up the change and create a hard link in the Hub's data archive directory to the
newly produced file. Then the uploader will take care of uploading it to the
MongoDB configured.

### Setting up a new source

Copy over an existing source, and modify the `Dockerfile` as you wish. Nothing
is strictly necessary besides setting up the correct user and permissions. But
it is recommended to keep `cron` installed, and reuse `docker-entrypoint.sh`.
Re-using the `run-spider.sh` is also recommended, but you don't need to use
Scrapy. Write your code to do the crawling and have it save its results to
`/data/<source_name>_crawled/data.ndjson` and save the UTC date and time of
at start of crawling, in ISO-8601 to `/data/<source_name>_crawled/release.txt`.
Again, do not  simply overwrite the files, write to temporary files and rename
them in the end. Build and test your code in a new container, and if that works,
include it in the `docker-compose.yml` file.

Then add a new source under `biothings-hub/files/nde-hub/hub/dataload/sources`,
create a Python package with `__init__.py`, `dumper.py`, and `uploader.py`. Do
the usual BioThings business there. I would recommend subclassing 
`biothings-hub/files/nde-hub/hub/dataload/nde.NDEFileSystemDumper` and 
`biothings-hub/files/nde-hub/hub/dataload/nde.NDESourceUploader` so you only
need to set `SRC_NAME` and `name` as class variables. Import your dumper and 
uploader in the newly created `__init__.py`. Re-build and re-start the
`biothings-hub` container with your new image. Everything should be ready.

## Current Issues
1. Due to the data transformation being incomplete. While data can be dumped and
uploaded in the Hub, it cannot be indexed to Elasticsearch as it conflicts with
the given mapping.
2. Cron does not forward signals and does not seem to respond to signals when 
running as PID 1 process. It is a minor issue but the crawler does not know when
it is about to get terminated, so it can't respond to it, like deleting the
temporary files. There may also be a very slight chance that it was in the
process of renaming the files so you get one file that is newer but the other is
not. The file renaming in the Scrapy NDJSON pipeline was designed with this flaw
in mind so the data file is renamed before the release information. Combined
with how the NDE Dumper only updates when the release is new, this should not be
an issue. However, in the long run, cruft temporary files may pile up in the
crawler output directories. Just delete them if they become a problem.
