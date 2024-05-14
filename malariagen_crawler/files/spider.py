# -*- coding: utf-8 -*-
# need to specify pythonpath because scrapy expects you to have the standard scrapy project structure
# need to specify settings because there is no scrapy.cfg file
# PYTHONPATH="$PYTHONPATH:." SCRAPY_SETTINGS_MODULE="settings" scrapy runspider spider.py
# If you are using the Dockerfile, runspider does this for you: /home/biothings/run-spider.sh

import logging
from datetime import datetime

import scrapy

logger = logging.getLogger("nde-logger")



# may take some time to start up as getting the ids takes a while
class MalariaGenSpider(scrapy.Spider):
    name = "malariagen"

    custom_settings = {
        "ITEM_PIPELINES": {
            "pipeline.MalariaGenItemProcessorPipeline": 100,
            "ndjson.NDJsonWriterPipeline": 999,
        }
    }

    start_urls = ['https://www.malariagen.net/data/archive/#']

    def parse(self, response):
        # Extract dataset links
        for a in response.css('div.card-title > a'):
            url = response.urljoin(a.attrib['href'])
            yield response.follow(url, self.parse_dataset)


    def parse_dataset(self, response):
        # Extract the release date
        release_text = response.xpath('//div[@class="hero-summary standard white"]/p[contains(text(), "Released on")]/text()').get()
        date_published = 'Not available'
        if release_text:
            date_published_raw = release_text.replace("Released on", "").strip().strip('.')
            try:
                date_published = datetime.strptime(date_published_raw, '%d %b %Y').strftime('%Y-%m-%d')
            except ValueError:
                pass


        # Distribution information
        distribution = []
        for block in response.css('blockquote.data_set-item-block-card'):
            title = block.css('h4.card-title::text').get().strip()
            description_texts = ' '.join(block.css('p::text').getall()).strip()
            href_url = response.urljoin(block.css('a::attr(href)').get(default=None))
            distribution.append({
                'name': title,
                'url': href_url,
                'description': description_texts
            })

        # Author information
        author_name = response.xpath('//blockquote[contains(.//h3, "Data package contact")]//a/text()').get(default=None)
        author_url = response.urljoin(response.xpath('//blockquote[contains(.//h3, "Data package contact")]//a/@href').get(default=None))

        # Initialize description collection
        description_parts = []

        # Use XPath to iterate through elements before the 'Data sets' header
        elements = response.xpath("//div[contains(@class, 'container container-2-col txt-block')]/*")
        for element in elements:
            # Check if we have reached the 'Data sets' heading
            if element.xpath('self::h3[contains(text(), "Data sets")]'):
                break  # Stop collecting text if 'Data sets' heading is reached

            # Otherwise, continue collecting text
            texts = element.xpath('.//text()').getall()
            cleaned_text = ' '.join(texts).replace('\n', ' ').strip()
            description_parts.append(cleaned_text)

        # Join all parts to form the full description
        description = ' '.join(description_parts).strip()

        # Extracting keywords dynamically from buttons
        bad_keywords = ["Download", "Go to"]
        keywords = response.css('a.button[class*="button-"][class*="-circle"]::text').getall()
        # remove bad keywords if they start with "Download PDF" or "Go to"
        keywords = [keyword for keyword in keywords if not any(keyword.startswith(bad_keyword) for bad_keyword in bad_keywords)]

        # Extracting the correct citation URL
        citation_url = response.xpath('//blockquote[contains(.//h3, "Citations")]//a[@data-wpel-link="external"]/@href').get(default=None)

        # Extracting 'isPartOf' data correctly, now capturing all text within the <a> tag, including nested tags
        project_name_parts = response.css('.hero-info.white a *::text').getall()  # Extracts all text components including from within <i>
        project_name = ''.join(project_name_parts).strip()  # Joins all parts together into one string
        project_url = response.urljoin(response.css('.hero-info.white a::attr(href)').get(default=None))

        open_urls = ["https://www.malariagen.net/project/pf3k/", "https://www.malariagen.net/project/p-falciparum-community-project/"]
        if project_url in open_urls:
            access_text = "Open"
        else:
            access_text = "Restricted"


        details = {
            'url': response.url,
            'name': response.css('title::text').get().strip(),
            'isPartOf': {
                '@type': "ResearchProject",
                'name': project_name,
                'url': project_url
            },
            'datePublished': date_published,
            'keywords': keywords,
            'conditionsOfAccess': access_text,
            'description': description,
            'distribution': distribution,
            'author': {
                'name': author_name,
                'url': author_url
            },
            'citation': response.urljoin(citation_url)
        }
        yield details
