# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: https://docs.scrapy.org/en/latest/topics/item-pipeline.html


# useful for handling different item types with a single interface
# mapping https://docs.google.com/spreadsheets/d/19hQf4sQ6ZLwvt8ADYyvooYrqmQCsxh386kYnycX3xmg/edit#gid=0

import datetime
import json
import logging
import re

logger = logging.getLogger("nde-logger")


__all__ = [
    "BeiItemProcessorPipeline",
]


def insert_value(d, key, value, extend=False):
    """ Insert a value into a dictionary, handling existing keys by converting to lists or extending strings as needed.
    """

    if key in d and not extend:
        if isinstance(d[key], list) and value not in d[key]:
            d[key].append(value)
        if not isinstance(d[key], list) and d[key] != value:
            d[key] = [d[key], value]
    elif d.get(key) and extend:
        d[key] = (d.get(key) + " " + value).strip()
    else:
        d[key] = value


class BeiItemProcessorPipeline:

    def process_item(self, item):
        _id = item.get("_id")
        url = item.get("url")
        output = {
            "@context": "http://schema.org/",
            "@type": "Dataset",
            "_id": "bei_" + _id.casefold(),
            "identifier": _id,
            "url": url,
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "BEI Resources",
                "url": "https://www.beiresources.org/Home.aspx",
                "versionDate": datetime.date.today().isoformat(),
                "archivedAt": url,
            },
        }

        if name := item.get("Organism:"):
            insert_value(output, "species", {"name": name})

        if "Product Name:" not in item and "Designations:" in item:
            name = item.get("Designations:")
            insert_value(output, "name", name)

        if keywords := item.get("Biosafety Level:"):
            try:
                keywords = int(keywords[0])
                insert_value(output, "keywords", f"Biosafety Level {keywords}")
            except (ValueError, TypeError):
                pass

        if sample_availability := item.get("Availability Status:"):
            if sample_availability.strip().casefold() == "in stock" or sample_availability.strip().casefold() == "made to order":
                insert_value(output, "sample", {"sampleAvailability": True})
            else:
                insert_value(output, "sample", {"sampleAvailability": False})

        if name := item.get("Storage Temperature:"):
            insert_value(output.setdefault("sample", {}), "sampleStorageTemperature", {"name": name})

        if name := item.get("Contributor:"):
            insert_value(output, "contributor", {"name": name})

        if description := item.get("Comments:"):
            insert_value(output, "description", description, extend=True)

        if credit_text := item.get("Citations:"):
            insert_value(output, "creditText", credit_text, extend=True)

        if citation := item.get("Publications Citing this Reagent:"):
            insert_value(output, "citedBy", {"citation": citation})

        if license := item.get("Restrictions On Use:"):
            insert_value(output, "license", license)

        if name := item.get("Product Name:"):
            insert_value(output, "name", name)

        if description := item.get("Description:"):
            insert_value(output.setdefault("usageInfo", {}), "description", description, extend=True)

        if credit_text := item.get("References:"):
            insert_value(output, "creditText", credit_text, extend=True)

        if description := item.get("Description:"):
            insert_value(output, "description", description, extend=True)

        if sample_state := item.get("Material Provided:"):
            insert_value(output.setdefault("sample", {}), "sampleState", sample_state, extend=True)

        if name := item.get("Manufacturer:"):
            insert_value(output, "author", {"@type": "Organization", "name": name})

        if description := item.get("Product Description:"):
            insert_value(output, "description", description, extend=True)

        if sample_process := item.get("Packing/Storage:"):
            insert_value(output.setdefault("sample", {}), "sampleProcess", sample_process)

        if experimental_purpose := item.get("Functional Activity:"):
            insert_value(output.setdefault("sample", {}), "experimentalPurpose", experimental_purpose)

        if citation := item.get("References:"):
            insert_value(output, "isBasedOn", {"citation": citation})

        if identifiers := item.get("Components:"):
            identifiers = identifiers.split(",")
            identifiers = [id.strip() for id in identifiers if id.strip()]
            output["hasPart"] = [{"identifier": id} for id in identifiers]

        if name := item.get("Additional Information:"):
            insert_value(output, "isRelatedTo", {"name": name})

        if description := item.get("Ownership statement:"):
            insert_value(output.setdefault("usageInfo", {}), "description", description, extend=True)

        if additional_property := item.get("Insert Size:"):
            insert_value(output.setdefault("sample", {}), "additionalProperty", {"@type": "PropertyValue", "name": "Insert Size", "value": additional_property})

        if not "Organism:" in item and "Taxonomy:" in item:
            name = item.get("Taxonomy:")
            insert_value(output, "species", {"name": name})

        if additional_property := item.get("Growth Conditions:"):
            insert_value(output.setdefault("sample", {}), "additionalProperty", {"@type": "PropertyValue", "name": "Growth Conditions", "value": additional_property})

        if additional_property := item.get("Safety Precautions:"):
            insert_value(output.setdefault("sample", {}), "additionalProperty", {"@type": "PropertyValue", "name": "Safety Precautions", "value": additional_property})

        if additional_property := item.get("Thawing and Growth:"):
            insert_value(output.setdefault("sample", {}), "additionalProperty", {"@type": "PropertyValue", "name": "Thawing and Growth", "value": additional_property})

        if description := item.get("Patents or other restrictions:"):
            insert_value(output.setdefault("usageInfo", {}), "description", description, extend=True)

        if additional_property := item.get("SubTaxon:"):
            insert_value(output.setdefault("species", {}), "additionalProperty", {"@type": "PropertyValue", "name": "SubTaxon", "value": additional_property})

        if description := item.get("Additional Notes:"):
            insert_value(output, "description", description, extend=True)

        if additional_property := item.get("Reconstitution:"):
            insert_value(output.setdefault("sample", {}), "additionalProperty", {"@type": "PropertyValue", "name": "Reconstitution", "value": additional_property})

        if sample_state := item.get("Solubility:"):
            insert_value(output.setdefault("sample", {}), "sampleState", sample_state, extend=True)

        if sample_state := item.get("Storage of Reconstituted Peptides:"):
            insert_value(output.setdefault("sample", {}), "sampleState", sample_state, extend=True)

        if additional_property := item.get("Sequence:"):
            insert_value(output.setdefault("sample", {}), "additionalProperty", {"@type": "PropertyValue", "name": "Sequence", "value": additional_property})

        return output

