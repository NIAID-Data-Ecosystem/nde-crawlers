import datetime
import os
import platform

import orjson
import scrapy.crawler
import scrapy.signals
import scrapy.spiders
from itemadapter import ItemAdapter


class NDJsonWriterPipeline:
    @classmethod
    def from_crawler(cls, crawler: scrapy.crawler.Crawler):
        # so that when spider_error occurs, data is not written
        ndjson_writer = cls()
        crawler.signals.connect(
            receiver=ndjson_writer.invalidate_data,
            signal=scrapy.signals.spider_error
        )
        return ndjson_writer

    def invalidate_data(self):
        self._valid = False

    def open_spider(self, spider: scrapy.spiders.Spider):
        # set the release string to be ISO date format
        # minute precision is good enough but feel free to change
        self.release_string = datetime.datetime.now(
            datetime.timezone.utc
        ).strftime('%Y-%m-%dT%H:%M:%SZ')
        dirname = os.path.join(
            '/data', f'{spider.name}_crawled'
        )
        os.makedirs(dirname, exist_ok=True)
        self.release_filename = os.path.join(
            dirname, 'release.txt'
        )
        # so that updates are as atomic as possible, using rename
        self.final_data_filename = os.path.join(
            dirname, 'data.ndjson'
        )
        self.tmp_filename = f'{self.final_data_filename}.{platform.node()}.{os.getpid()}'
        self.rl_tmp_filename = f'{self.release_filename}.{platform.node()}.{os.getpid()}'
        with open(self.rl_tmp_filename, 'w') as release_file:
            release_file.write(self.release_string)
        os.makedirs(os.path.dirname(self.tmp_filename), exist_ok=True)
        self.fd = open(self.tmp_filename, 'wb')
        self._valid = True

    def close_spider(self, spider: scrapy.spiders.Spider):
        self.fd.close()
        # data is invalid, just remove and leave the scene
        if not self._valid:
            spider.logger.warning(
                "Errors occurred while running, so not saving potentially corrupt data."
            )
            os.unlink(self.tmp_filename)
            return
        try:
            # it is critical that we use rename here so that
            # the files have a different inode
            # because NDE Dumper uses hardlinks
            os.rename(self.tmp_filename, self.final_data_filename)
            os.rename(self.rl_tmp_filename, self.release_filename)
        except:
            # cleanup the mess by attempting to remove everything
            # so that bad data does not propagate
            spider.logger.error(
                "Error updating data/release, will remove corrupt files")
            errors = []
            for fn in [
                self.tmp_filename,
                self.final_data_filename,
                self.rl_tmp_filename,
                self.release_filename
            ]:
                try:
                    os.unlink(fn)
                except FileNotFoundError:
                    pass
                except Exception as e:
                    spider.logger.critical(
                        f"Failed to remove corrupt file {fn}")
                    errors.append(e)
            if errors:
                # re-raise all errors
                raise RuntimeError(errors)

    def process_item(self, item, spider):
        line = orjson.dumps(ItemAdapter(item).asdict()) + b"\n"
        self.fd.write(line)
        return item
