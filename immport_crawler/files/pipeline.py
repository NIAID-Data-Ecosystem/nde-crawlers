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
        if _id := item.pop('_id', None):
            item['_id'] = 'IMMPORT_' + _id.rsplit('/', 1)[-1]
        if author := item.pop('creator', None):
            for data in author:
                for key,value in data.items():
                    data[key] = {"name": value} if key == "affiliation" else value
            item['author'] = author
        if cited_by := item.pop('citations', None):
            item['citedBy'] = cited_by
        if identifier := item.pop('identifiers', None):
            item['identifier'] = identifier
        if species := item.pop('species', None):
            item['species'] = {'name': species}
        if measurement_technique := item.pop('measurementTechnique', None):
            item['measurementTechnique'] = {'name': measurement_technique}
        if health_conditions := item.pop('keywords', None):
            no_dups = set()
            hc = []
            for health_condition in health_conditions:
                lower = health_condition.lower()
                if lower not in no_dups:
                    no_dups.add(lower)
                    hc.append({'name': lower})
            item['healthCondition'] = hc

        if date := item.pop('curationDate', None):
            date = datetime.datetime.strptime(date, "%m/%d/%Y")  # mm/dd/YYYY is my guess
            if distribution := item.pop('distribution', None):
                item['distribution'] = _set_single_element_or_all_in_list(
                    distribution, 'dateModified', date
                )
            if curated_by := item.pop('curatedBy', None):
                # remove curatedBy
                # item['curatedBy'] = _set_single_element_or_all_in_list(
                #     curated_by, 'curationDate', date
                # )
                pass
            if item.get('includedInDataCatalog'):
                item['includedInDataCatalog']['versionDate'] = date
                item['includedInDataCatalog']['name'] = "ImmPort"
            if 'date' not in item:
                item['date'] = date
        
        # temp fix for outbreak.info
        if cb_outbreak := item.get('includedInDataCatalog'):
            item['curatedBy'] = cb_outbreak
        return item
