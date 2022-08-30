# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: https://docs.scrapy.org/en/latest/topics/item-pipeline.html


# useful for handling different item types with a single interface

import datetime
import logging
import re

logger = logging.getLogger('nde-logger')


__all__ = [
    'DryadItemProcessorPipeline',
]


class DryadItemProcessorPipeline:

    def format_date(self, date_string):
        """ Formats python and military datetime into an isoformat date

        Args:
            da  te_string: a date string
        Returns: An isoformatted date if there is a datetime if not then return None
        """
        if re.match("\d+-\d+-\d+", date_string):
            if 'T' in date_string:
                date = datetime.datetime.fromisoformat(date_string.split('T')[0]).date().isoformat()
            else:
                date = datetime.datetime.fromisoformat(date_string).date().isoformat()
            return date
        return None

    def process_item(self, item: dict, spider):
        output = {
            "@context": item.pop('@context'),
            "@type": item.pop('@type'),
            "url": item.pop('url'),
            "_id": "DRYAD_" + item.pop('@id').split('//')[-1].replace('/', ':'),
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "Dryad Digital Repository",
                "url": "https://datadryad.org",
                'versionDate': datetime.date.today().isoformat()
            }
        }

        if name := item.pop('name', None):
            output['name'] = name

        if description := item.pop('description', None):
            output['description'] = description

        if content_url := item.pop('contentUrl', None):
            output['contentUrl'] = content_url

        if identifier := item.pop('identifier', None):
            output['identifier'] = identifier
            if "doi" in identifier:
                output['doi'] = identifier

        if keywords := item.pop('keywords', None):
            output['keywords'] = keywords

        if author := item.pop('creator', None):
            output['author'] = author

        if distribution := item.pop('distribution', None):
            output['distribution'] = distribution

        # There are 2 different cases and 4 different time format
        """
            Case 1: No list (contains one of the formats)
            Case 2: List (we assume it contains format 1 and  at least one format 2-4)
            Format 1: Y (2019)
            Format 2: Y-M-D (2018-12-13)
            Format 3: Python datetime (2015-10-16T00:00:00+00:00)
            Format 4: Military datetime (2015-10-16T00:00:00Z)
        """
        if temp_cov := item.pop('temporalCoverage', None):
            if isinstance(temp_cov, list):
                for tc in temp_cov:
                    # we assume there will be a format 2-4 in a list
                    if date := self.format_date(tc):
                        output['temporalCoverage'] = {'temporalInterval': {'endDate': date}}
            else:
                # we cannot assume format 2-4 is in the no list
                if date := self.format_date(temp_cov):
                    pass
                else:
                    date = datetime.datetime.strptime(temp_cov, '%Y').date().isoformat()
                output['temporalCoverage'] = {'temporalInterval': {'endDate': date}}

            if not output.get('temporalCoverage'):
                logger.info('No temporal coverage: ' + " ".join(temp_cov))
                logger.info('ID is %s', output['_id'])

        if spatial_covs := item.pop('spatialCoverage', None):
            if isinstance(spatial_covs, list):
                sc = []
                for spatial_cov in spatial_covs:
                    if isinstance(spatial_cov, str):
                        sc.append({'name': spatial_cov})
                output['spatialCoverage'] = sc
            else:
                if isinstance(spatial_covs, str):
                    output['spatialCoverage'] = {'name': spatial_covs}

        if citation := item.pop('citation', None):
            output['citation'] = {'url': citation}

        if license_obj := item.pop('license', None):
            output['license'] = license_obj['license']

        if date_pub := item.pop('datePublished'):
            date_pub = datetime.datetime.strptime(date_pub.split(': ')[1], '%B %d, %Y').date().isoformat()
            output['datePublished'] = date_pub

        item.pop('provider', None)
        item.pop('isAccessibleForFree', None)
        item.pop('publisher', None)
        item.pop('version', None)

        if item:
            logger.warning("Haven't parsed all keys in dryad_crawler: " + "\t".join(item.keys()))
        
        return output
