import datetime


__all__ = [
    'Immport2OutbreakDatasetPipeline',
]


def _set_single_element_or_all_in_list(obj, key, value):
    """mutates obj, list of dict or dict"""
    if isinstance(obj, list):
        for element in obj:
            element[key] = value
    else:
        obj[key] = value
    return obj


class Immport2OutbreakDatasetPipeline:
    def process_item(self, item: dict, _):
        # FIXME: the format is surely STILL NOT VALID
        #  but I'll leave it here for someone to fix the transformation
        # see discussion
        #  https://suwulab.slack.com/archives/C01078YR6SW/p1631035614003000?thread_ts=1630687426.004000&cid=C01078YR6SW
        if author := item.pop('creator', None):
            item['author'] = author
        if cited_by := item.pop('citation', None):
            item['citedBy'] = cited_by
        if date := item.pop('curationDate', None):
            date = datetime.datetime.strptime(date, "%m/%d/%Y")  # mm/dd/YYYY is my guess
            if distribution := item.pop('distribution', None):
                item['distribution'] = _set_single_element_or_all_in_list(
                    distribution, 'dateModified', date
                )
            if curated_by := item.pop('curatedBy', None):
                item['curatedBy'] = _set_single_element_or_all_in_list(
                    curated_by, 'curationDate', date
                )
            if 'date' not in item:
                item['date'] = date
        return item
