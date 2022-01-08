import datetime
import urllib.parse
import re

import orjson
import requests
import scrapy.core.scraper
from scrapy.http import Request, Response


class ZenodoSpider(scrapy.Spider):
    """
    Crawl Zenodo

    Loops through their API, accounting for the throttling

    Zenodo API can only read 10k records from single query so we group
    by publication date

    we bump the range up or down so we aim for approx. 5k - 10k documents
    per query
    """

    name = 'zenodo'
    custom_settings = {
        'ITEM_PIPELINES': {
            'pipeline.ZenodoItemProcessorPipeline': 100,
            'ndjson.NDJsonWriterPipeline': 999,
        }
    }
    max_docs = 10000

    def start_requests(self):
        self.default_params = {
            'page': '1',
            'size': '1000',
        }
        self.base_url = 'https://zenodo.org/api/records/'
        self.date_interval = 365 * 10  # start with 10 years, we adjust as needed
        # find the last publication date
        params = {
            'sort': '-publication_date',
            'size': '1',
            'page': '1',
        }
        fd_resp = requests.get(self.base_url, params=params).json()
        self.final_date = datetime.date.fromisoformat(
            fd_resp['hits']['hits'][0]['metadata']['publication_date']
        )
        # for debugging, uncomment to only scrape a few docs.
        # self.final_date = datetime.date.fromisoformat('1800-01-01')
        self.logger.info("last publication date is %s", self.final_date)

        # we are betting that almost no documents hit so this query is good enough
        params = {
            'q': '_missing_:publication_date',
        }
        params.update(self.default_params)
        yield Request(
            f'{self.base_url}?' + urllib.parse.urlencode(params),
            callback=self.parse_no_date
        )
        start_date = datetime.date.fromisoformat('0001-01-01')
        end_date = datetime.date.fromisoformat('1700-01-01')
        # end_date is heuristics, change if needed or figure out something better
        yield self._produce_date_req(start_date, end_date)

    def _parse_body(self, resp: dict, callback=None):
        for hit in resp['hits']['hits']:
            yield hit
        if next_page := resp.get('links', {}).get('next', None):
            yield Request(next_page, callback=callback)

    def parse_no_date(self, response: Response, **kwargs):
        # parse where query does not involve dates, i.e. q=_missing_:publication_date
        resp = orjson.loads(response.body)
        # we are assuming less than 10k
        total = resp['hits']['total']
        if total > self.max_docs:
            self.logger.warning(
                "%s produced %d documents, can only fetch %d max.",
                response.url,
                total,
                self.max_docs
            )
        yield from self._parse_body(resp, callback=self.parse_no_date)

    def _produce_date_req(self, start_date: datetime.date, end_date: datetime.date):
        if start_date > self.final_date:
            return None  # this would be the ending point
        self.logger.debug(
            "Next query range is publication dates %s - %s inclusive",
            start_date,
            end_date,
        )
        params = {
            'q': f'publication_date: '
                 f'[{start_date.isoformat()} TO {end_date.isoformat()}]',
        }
        params.update(self.default_params)
        return Request(
            f'{self.base_url}?' + urllib.parse.urlencode(params)
        )

    def _calc_dates(self, curr_end_date):
        # calculate new dates, inclusive
        # so when self.date_interval == 1, start = end
        start_date = curr_end_date + datetime.timedelta(days=1)
        end_date = start_date + datetime.timedelta(days=self.date_interval - 1)
        return start_date, end_date

    def parse(self, response: Response, **kwargs):
        # this handles requests with dates in the query
        qs = urllib.parse.urlparse(response.url).query
        # figure out the range of the current query
        parsed_qs = urllib.parse.parse_qs(qs)
        dates = re.findall(
                    r'\d{4}-\d{2}-\d{2}',
                    parsed_qs['q'][0]
                )
        start_date = datetime.date.fromisoformat(dates[0])
        end_date = datetime.date.fromisoformat(dates[1])

        resp = orjson.loads(response.body)
        total = resp['hits']['total']
        if total <= 10000:
            yield from self._parse_body(resp)
            if 'next' not in resp['links']:
                 # move on to next date range
                if resp['hits']['total'] < 5000:
                    self.date_interval *= 2
                    self.logger.info(
                        "Bumping query date interval to %d days",
                        self.date_interval
                    )
                # calculate new dates
                start_date, end_date = self._calc_dates(end_date)
                if next_req := self._produce_date_req(start_date, end_date):
                    yield next_req

        else:  # got too many hits
            curr_page = int(parsed_qs['page'][0])
            if self.date_interval == 1:
                self.logger.warning("%s has too many publications", start_date)
            if curr_page > 1 or self.date_interval == 1:
                yield from self._parse_body(resp)
            else:
                self.date_interval = self.date_interval >> 1  # half
                self.logger.info(
                    "Dropping query date interval to %d days",
                    self.date_interval
                )
                start_date, end_date = self._calc_dates(
                    start_date - datetime.timedelta(days=1)
                )
                yield self._produce_date_req(start_date, end_date)
