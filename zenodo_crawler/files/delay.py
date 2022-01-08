import datetime
import logging

import scrapy
import scrapy.signals
from scrapy.crawler import Crawler
from scrapy.spiders import Spider
from scrapy.http import Request, Response

logger = logging.getLogger(__name__)


class XRateLimitDelay:
    """
    Throttle based on X-RateLimit header from Zenodo

    Copy pasted from the Auto Throttle plugin from Scrapy
    """
    def __init__(self, crawler: Crawler):
        self.crawler = crawler
        self.min_delay = crawler.settings.getfloat('DOWNLOAD_DELAY')
        self.reserve_rate = crawler.settings.getint('RATELIMIT_RESERVE', 5)
        logger.info(
            "XRateLimitDelay enabled with min delay of %f and reserve rate limit of %d",
            self.min_delay,
            self.reserve_rate,
        )
        crawler.signals.connect(self._spider_opened, signal=scrapy.signals.spider_opened)
        crawler.signals.connect(self._response_downloaded, signal=scrapy.signals.response_downloaded)

    @classmethod
    def from_crawler(cls, crawler: Crawler):
        return cls(crawler)

    def _spider_opened(self, spider: Spider):
        spider.download_delay = self.min_delay
        logger.debug(
            "spider %s opened, setting delay to %f",
            Spider,
            spider.download_delay
        )

    def _response_downloaded(self, response, request, spider):
        key, slot = self._get_slot(request, spider)
        latency = request.meta.get('download_latency')
        if latency is None or slot is None:
            return

        olddelay = slot.delay
        self._adjust_delay(slot, latency, response)
        diff = slot.delay - olddelay
        size = len(response.body)
        conc = len(slot.transferring)
        logger.debug(
            "slot: %(slot)s | conc:%(concurrency)2d | "
            "delay:%(delay)5d ms (%(delaydiff)+d) | "
            "latency:%(latency)5d ms | size:%(size)6d bytes",
            {
                'slot': key, 'concurrency': conc,
                'delay': slot.delay * 1000, 'delaydiff': diff * 1000,
                'latency': latency * 1000, 'size': size
            },
            extra={'spider': spider}
        )

    def _get_slot(self, request: Request, spider: Spider):
        key = request.meta.get('download_slot')
        return key, self.crawler.engine.downloader.slots.get(key)

    def _adjust_delay(self, slot, latency, response: Response):
        """Define delay adjustment policy"""

        remaining_count = response.headers.get('X-RateLimit-Remaining', 60)  # assume 60
        remaining_count = int(remaining_count)
        # not using time.time because behavior is OS dependent
        now = datetime.datetime.now(datetime.timezone.utc).timestamp()
        reset = response.headers.get('X-RateLimit-Reset', None)
        if reset is None:
            reset = now + 60
        else:
            reset = float(reset)
        if remaining_count <= self.reserve_rate:
            new_delay = reset - now
        else:
            new_delay = self.min_delay

        new_delay = max(new_delay, self.min_delay)

        slot.delay = new_delay
