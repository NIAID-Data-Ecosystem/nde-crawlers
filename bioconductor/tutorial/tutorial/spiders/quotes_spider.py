import scrapy


class QuotesSpider(scrapy.Spider):
    name = "quotes"

    def start_requests(self):
        urls = [
            'https://bioconductor.org/packages/3.15/bioc/html/a4Core.html',
            'https://bioconductor.org/packages/3.15/bioc/html/diffloop.html',
        ]
        for url in urls:
            yield scrapy.Request(url=url, callback=self.parse)

    def parse(self, response):
        filename = response.url.split("/")[-1]
        with open(filename, 'wb') as f:
            f.write(response.body)
        self.log(f'Saved file {filename}')
