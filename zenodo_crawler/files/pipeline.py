import datetime

__all__ = [
    'ZenodoItemProcessorPipeline',
]


def _set_single_element_or_all_in_list(obj, key, value):
    """mutates obj, list of dict or dict"""
    if isinstance(obj, list):
        for element in obj:
            element[key] = value
    else:
        obj[key] = value
    return obj


class ZenodoItemProcessorPipeline:
    def process_item(self, item: dict, _):
        # FIXME: this is very incomplete and for demonstration only
        output = {
            '@context': "https://schema.org/",
            'curatedBy': {
                '@type': 'Organization',
                'name': 'Zenodo',
                'url': 'https://zenodo.org/',
                'versionDate': datetime.date.today().isoformat()
            },
            '@id': 'https://doi.org/' + item['doi'],
            '@type': item['metadata']['resource_type']['type'],
            'name': item['metadata']['title'],
            'author': [],
            'description': item['metadata']['description'],
            'identifier': item['doi'],
            'dateModified': datetime.datetime.fromisoformat(
                item['updated']
            ).date().isoformat(),
        }
        # if 'conceptrecid' in item:
        #     output['_id'] = 'zenodo.' + item['conceptrecid']
        output['_id'] = 'zenodo.' + str(item['id'])
        for person in item['metadata']['creators']:
            if 'name' in person:
                output['author'].append(
                    {
                        '@type': 'Person',
                        'name': person['name'],
                    }
                )
        if 'publication_date' in item['metadata']:
            output['datePublished'] = item['metadata']['publication_date']
        return output
