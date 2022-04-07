import datetime
from biothings import config


# helper function to get the date transformation
def add_date(doc):
    dates = []
    if doc.get('date'):
        dates.append(doc.get('date'))
    if doc.get('dateCreated'):
        dates.append(doc.get('dateCreated'))
    if doc.get('dateModified'):
        dates.append(doc.get('dateModified'))
    if doc.get('datePublished'):
        dates.append(doc.get('datePublished'))
    if dates:
        dates.sort()
        date = datetime.datetime.fromisoformat(dates[-1]).date().isoformat()
        doc['date'] = date

    return doc
