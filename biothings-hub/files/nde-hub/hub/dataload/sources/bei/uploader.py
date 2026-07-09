import re

from hub.dataload.nde import NDESourceSampleUploader
from utils.extract import process_descriptions
from utils.pubtator import standardize_data
from utils.utils import nde_upload_wrapper


AMBIGUOUS_TAXON_RE = re.compile(r"\b(?:unknown|unidentified)\b", re.IGNORECASE)
TAXON_LABEL_FIELDS = ("name", "commonName", "displayName", "originalName")
TAXON_FIELDS = ("species", "infectiousAgent")


def _iter_taxon_labels(entry):
    if isinstance(entry, str):
        yield entry
        return

    if not isinstance(entry, dict):
        return

    for field in TAXON_LABEL_FIELDS:
        value = entry.get(field)
        if isinstance(value, str):
            yield value

    alternate_names = entry.get("alternateName")
    if isinstance(alternate_names, str):
        yield alternate_names
    elif isinstance(alternate_names, list):
        for alternate_name in alternate_names:
            if isinstance(alternate_name, str):
                yield alternate_name


def _is_ambiguous_taxon(entry):
    return any(AMBIGUOUS_TAXON_RE.search(label) for label in _iter_taxon_labels(entry))


def _filter_taxon_entries(value):
    entries = value if isinstance(value, list) else [value]
    filtered = [entry for entry in entries if not _is_ambiguous_taxon(entry)]
    if not filtered:
        return None
    return filtered if isinstance(value, list) else filtered[0]


def remove_ambiguous_taxonomy(doc):
    """Drop BEI sample taxonomy terms such as unknown or unidentified taxa."""
    if doc.get("@type") != "Sample":
        return doc

    for field in TAXON_FIELDS:
        if field not in doc:
            continue
        filtered = _filter_taxon_entries(doc[field])
        if filtered is None:
            doc.pop(field, None)
        else:
            doc[field] = filtered

    return doc


class BeiUploader(NDESourceSampleUploader):
    name = "bei"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = standardize_data(data_folder)
        docs = process_descriptions(docs)
        for doc in docs:
            remove_ambiguous_taxonomy(doc)
            yield doc
