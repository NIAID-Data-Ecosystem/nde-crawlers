# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: https://docs.scrapy.org/en/latest/topics/item-pipeline.html

import datetime
import logging
import re

logger = logging.getLogger("nde-logger")


__all__ = [
    "CeirrItemProcessorPipeline",
]


INVALID_VALUES = {
    "",
    "-",
    "--",
    "na",
    "n/a",
    "none",
    "null",
    "not applicable",
    "not available",
    "unknown",
    "unspecified",
}

SEGMENT_FIELDS = ["HA", "NA", "NP", "NS", "PA", "PB1", "PB2", "MP"]

CEIRR_CATALOG = {
    "@type": "DataCatalog",
    "name": "Centers of Excellence for Influenza Research and Response (CEIRR) Resources",
    "alternateName": ["CEIRR Resources", "CEIRR Network Resources"],
    "identifier": "dde_5a908221f16d72c6",
    "url": "https://www.ceirr-network.org/",
}

CEIRR_USAGE_INFO = {
    "@type": "CreativeWork",
    "description": (
        "The CEIRR Reagents webpage is a searchable and filterable catalog of over 10,000 free "
        "and unique reagents developed by CEIRR investigators, which are made available to the "
        "broad scientific community to help advance influenza and other viral research worldwide. "
        "Some reagents on this page may require additional documentation to verify that the "
        "requesting lab has the appropriate permits and facilities to receive the reagent. "
        "Additional time may be required to ship reagents that are classified as BSL-3. Please "
        "respond promptly to any correspondence to ensure that the reagent is shipped as quickly "
        "as possible."
    ),
}

NIAID_FUNDER = {
    "@type": "Organization",
    "name": "National Institute of Allergy and Infectious Diseases",
    "alternateName": "NIAID",
    "url": "https://www.niaid.nih.gov",
}

CEIRR_FUNDING_IDS = [
    "75N93021C00007",
    "75N93021C00014",
    "75N93021C00015",
    "75N93021C00017",
    "75N93021C00045",
    "OD037684",
    "75N93021C00018",
    "75N93021C00016",
]

CEIRR_SOURCE_ORGANIZATION = {
    "@type": "ResearchProject",
    "name": "NIAID CEIRR Network",
    "abstract": (
        "The Centers of Excellence for Influenza Research and Response (CEIRR) is an "
        "international research network that studies the transmission and pathogenesis of "
        "influenza while providing critical infrastructure for outbreak responses. The program "
        "builds on prior influenza research initiatives to enhance global preparedness strategies."
    ),
    "description": (
        "NIAID established the Centers of Excellence for Influenza Research and Response (CEIRR) "
        "to study natural history, transmission, and pathogenesis of influenza and provide an "
        "international research infrastructure to address influenza outbreaks. CEIRR replaced the "
        "Centers of Excellence for Influenza Research and Surveillance (CEIRS) program, which was "
        "supported by contracts that concluded on March 31, 2021. For more information, visit the "
        "NIAID program page: https://www.niaid.nih.gov/research/centers-excellence-influenza-research-response"
    ),
    "alternateName": ["CEIRR", "Centers of Excellence for Influenza Research and Response"],
    "url": "https://www.ceirr-network.org/",
    "parentOrganization": ["NIAID"],
}

TOPIC_CATEGORIES = [
    {
        "@type": "DefinedTerm",
        "name": "Sample collections",
        "identifier": "topic_3277",
        "url": "http://edamontology.org/topic_3277",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Infectious disease",
        "identifier": "topic_3324",
        "url": "http://edamontology.org/topic_3324",
        "inDefinedTermSet": "EDAM",
    },
    {
        "@type": "DefinedTerm",
        "name": "Immunology",
        "identifier": "topic_0804",
        "url": "http://edamontology.org/topic_0804",
        "inDefinedTermSet": "EDAM",
    },
]

PROJECT_NAMES = {
    "SJCEIRR": "St Jude Center of Excellence for Influenza Research and Response",
    "SJCEIRS": "St Jude Centers of Excellence for Influenza Research and Surveillance",
    "Penn-CEIRR": "University of Pennsylvania Center of Excellence for Influenza Research and Response",
    "JHCEIRR": "Johns Hopkins Center of Excellence for Influenza Research and Response",
    "Emory-CEIRR": "Emory University Center of Excellence for Influenza Research and Response",
    "CRIPT": "Center for Research on Influenza Pathogenesis and Transmission",
    "CIDER": "Center for Influenza Disease and Emergence Research",
}


def insert_value(d, key, value, extend=False):
    """Insert a value, preserving distinct values as lists."""

    if value is None:
        return

    if key in d and not extend:
        if isinstance(d[key], list):
            if value not in d[key]:
                d[key].append(value)
        elif d[key] != value:
            d[key] = [d[key], value]
    elif d.get(key) and extend:
        d[key] = (d.get(key) + " " + value).strip()
    else:
        d[key] = value


def clean_value(value):
    if value is None:
        return None
    value = str(value).replace("\xa0", " ")
    value = re.sub(r"\s+", " ", value).strip()
    if value.casefold() in INVALID_VALUES:
        return None
    return value


def maybe_number(value):
    value = clean_value(value)
    if value is None:
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    return int(number) if number.is_integer() else number


def date_only(value):
    value = clean_value(value)
    if not value:
        return None
    # CEIRR emits ISO-like timestamps. Keeping the date avoids adding a parser
    # dependency while still fitting the hub date mapping.
    if re.match(r"^\d{4}-\d{2}-\d{2}", value):
        return value[:10]
    return None


def year_only(value):
    value = clean_value(value)
    if value and re.fullmatch(r"\d{4}", value):
        return value
    return None


def as_identifier_list(values):
    identifiers = []
    seen = set()
    for value in values:
        value = clean_value(value)
        if not value:
            continue
        key = value.casefold()
        if key not in seen:
            identifiers.append(value)
            seen.add(key)
    if not identifiers:
        return None
    return identifiers[0] if len(identifiers) == 1 else identifiers


def add_additional_property(output, property_id, value, name=None):
    value = clean_value(value)
    if not value:
        return
    property_value = {"@type": "PropertyValue", "propertyID": property_id, "value": value}
    if name:
        property_value["name"] = name
    insert_value(output, "additionalProperty", property_value)


def expand_abbreviation(item, group, value):
    value = clean_value(value)
    if not value:
        return None

    abbreviations = item.get("_ceirr_abbreviations") or {}
    group_values = abbreviations.get(group) or {}
    if value in group_values:
        return clean_value(group_values[value])

    if "-" in value:
        prefix, suffix = value.split("-", 1)
        if prefix in group_values:
            return clean_value(suffix) or clean_value(group_values[prefix])

    return value


def split_pmids(value):
    value = clean_value(value)
    if not value:
        return []
    return [pmid for pmid in re.split(r"[,;]\s*|\s+", value) if re.fullmatch(r"\d+", pmid)]


def split_multi_value(value):
    value = clean_value(value)
    if not value:
        return []
    values = [part.strip() for part in re.split(r"\s*;\s*|\s*,\s*", value) if part.strip()]
    return values or [value]


def normalize_unit(unit):
    unit = clean_value(unit)
    if not unit:
        return None
    normalized = unit.replace("\u03bc", "u").replace("\u00b5", "u").casefold()
    return {
        "ml": "mL",
        "ul": "uL",
        "l": "L",
    }.get(normalized, unit)


def add_sample_quantity(output, name, value=None, unit_text=None):
    quantity = {"name": name}
    if value is not None:
        quantity["value"] = value
    if unit_text:
        quantity["unitText"] = unit_text
    insert_value(output, "sampleQuantity", quantity)


def add_amount_quantity(output, text):
    text = clean_value(text)
    if not text:
        return

    match = re.search(r"(?P<value>\d+(?:\.\d+)?)\s*(?P<unit>mL|ml|ML|uL|ul|UL|\u00b5L|\u03bcL|L)\b", text)
    if match:
        value = maybe_number(match.group("value"))
        add_sample_quantity(output, text, value=value, unit_text=normalize_unit(match.group("unit")))
    else:
        add_sample_quantity(output, text)


def add_concentration_quantity(output, text):
    text = clean_value(text)
    if not text:
        return

    value = None
    unit_text = None
    sci_match = re.search(r"(?P<base>\d+(?:\.\d+)?)\s*[xX]\s*10(?:\^|E)?(?P<exp>[+-]?\d+)", text)
    if sci_match:
        value = float(sci_match.group("base")) * (10 ** int(sci_match.group("exp")))
    else:
        numeric_match = re.search(r"\d+(?:\.\d+)?(?:[eE][+-]?\d+)?", text)
        if numeric_match:
            value = maybe_number(numeric_match.group(0))

    if re.search(r"/\s*(mL|ml)\b", text):
        unit_text = "units/mL"

    add_sample_quantity(output, text, value=value, unit_text=unit_text)
    add_additional_property(output, "Concentration", text)


class CeirrItemProcessorPipeline:
    SOURCE_NAME = CEIRR_CATALOG["name"]
    SOURCE_URL = "https://www.ceirr-network.org/resources/reagents"
    BEI_HOME_URL = "https://www.beiresources.org/Home.aspx"

    def process_item(self, item):
        catalog_number = clean_value(item.get("_uuid")) or clean_value(item.get("sample_id"))
        if not catalog_number:
            raise ValueError(f"CEIRR record without a catalog number: {item}")

        source_url = clean_value(item.get("_ceirr_source_url")) or self.SOURCE_URL
        identifiers = [
            catalog_number,
            item.get("sample_identifier"),
            item.get("sample_id"),
            item.get("SK"),
        ]

        bei_number = clean_value(item.get("bei_number"))
        if bei_number:
            identifiers.extend([bei_number, self.bei_overlap_id(bei_number)])

        catalog = dict(CEIRR_CATALOG)
        catalog["versionDate"] = datetime.date.today().isoformat()
        catalog["archivedAt"] = source_url

        output = {
            "@context": "http://schema.org/",
            "@type": "Sample",
            "_id": "ceirr_" + catalog_number.casefold(),
            "identifier": as_identifier_list(identifiers),
            "url": source_url,
            "includedInDataCatalog": catalog,
            "conditionsOfAccess": "Restricted",
            "isAccessibleForFree": True,
            "additionalType": "BioSample",
            "usageInfo": dict(CEIRR_USAGE_INFO),
            "funding": [
                {"@type": "MonetaryGrant", "identifier": funding_id, "funder": dict(NIAID_FUNDER)}
                for funding_id in CEIRR_FUNDING_IDS
            ],
            "sourceOrganization": [dict(CEIRR_SOURCE_ORGANIZATION)],
            "topicCategory": [dict(topic) for topic in TOPIC_CATEGORIES],
        }

        if name := clean_value(item.get("reagent_name")) or clean_value(item.get("sample_identifier")):
            insert_value(output, "name", name)

        if alternate_id := clean_value(item.get("sample_identifier")):
            insert_value(output, "alternateIdentifier", alternate_id)

        if date_created := date_only(item.get("_sample_create_date")) or date_only(item.get("_create_timestamp")):
            output["dateCreated"] = date_created

        if date_modified := date_only(item.get("_modified_timestamp")) or date_only(item.get("modified_timestamp")):
            output["dateModified"] = date_modified

        if availability := clean_value(item.get("availability")):
            if availability.casefold() == "y":
                output["creativeWorkStatus"] = "Available"
                output["sampleAvailability"] = True
            else:
                output["sampleAvailability"] = False

        if category := clean_value(item.get("_ceirr_reagent_category")):
            insert_value(output, "sampleType", {"name": category})
            insert_value(output, "keywords", category)

        self.add_project(output, item)
        self.add_contributors(output, item)
        self.add_species(output, item)
        self.add_infectious_agent(output, item)
        self.add_biological_context(output, item)
        self.add_publications(output, item)
        self.add_bei_link(output, bei_number)
        self.add_gene_segments(output, item)

        return output

    def bei_overlap_id(self, bei_number):
        return f"bei_{bei_number.casefold()}"

    def add_project(self, output, item):
        project_id = clean_value(item.get("project_identifier"))
        if project_id:
            project_name = PROJECT_NAMES.get(project_id, project_id)
            organization = {"@type": "Organization", "identifier": project_id, "name": project_name}
            insert_value(output, "collector", organization)
            insert_value(output, "author", organization)

        if study_id := clean_value(item.get("study_identifier")):
            insert_value(output, "isPartOf", {"@type": "CreativeWork", "identifier": study_id})

    def add_contributors(self, output, item):
        contact_name = clean_value(item.get("contact_name"))
        contact_email = clean_value(item.get("contact_email"))
        institution = clean_value(item.get("contributing_institution"))
        if contact_name or contact_email:
            contributor = {"@type": "Person"}
            if contact_name:
                contributor["name"] = contact_name
            if contact_email:
                contributor["email"] = contact_email
            if institution:
                contributor["affiliation"] = institution
            insert_value(output, "contributor", contributor)

        if institution:
            insert_value(output, "itemLocation", {"@type": "Place", "identifier": institution})

    def add_species(self, output, item):
        host = (
            clean_value(item.get("host_common_name"))
            or clean_value(item.get("_common_host"))
            or clean_value(item.get("host"))
        )
        origin = clean_value(item.get("origin"))
        if host and host != origin:
            species = {"name": host}
            if host_id := clean_value(item.get("host_identifier")):
                species["identifier"] = host_id
            insert_value(output, "species", species)

        if sex := clean_value(item.get("host_sex")):
            sex_map = {"M": "Male", "F": "Female"}
            if sex.casefold() != "u":
                output["sex"] = sex_map.get(sex, sex)

        if strain := clean_value(item.get("host_strain")):
            insert_value(output, "associatedGenotype", strain)

    def add_infectious_agent(self, output, item):
        for strain_name in split_multi_value(item.get("strain_name")):
            insert_value(output, "infectiousAgent", {"name": strain_name})
            insert_value(output, "associatedGenotype", strain_name)

        for parent_strain in split_multi_value(item.get("parent_strain_name")):
            insert_value(output, "infectiousAgent", {"name": parent_strain})
            insert_value(output, "associatedGenotype", parent_strain)

        influenza_type = clean_value(item.get("influenza_type"))
        subtype = clean_value(item.get("subtype"))
        agent_terms = [term for term in [influenza_type, subtype] if term]
        if agent_terms:
            insert_value(output, "infectiousAgent", {"name": "Influenza " + " ".join(agent_terms)})

        if influenza_type:
            add_additional_property(output, "Influenza Type", influenza_type)

        if origin := clean_value(item.get("origin")):
            insert_value(output, "locationOfOrigin", {"name": origin})

        if year := year_only(item.get("year")):
            output["dateCollected"] = year

    def add_biological_context(self, output, item):
        if material := expand_abbreviation(item, "sample_material", item.get("sample_material")):
            insert_value(output, "sampleType", {"name": material})
            insert_value(output, "cellType", {"name": material})

        if material_form := clean_value(item.get("sample_material_form")):
            insert_value(output, "sampleState", material_form, extend=True)

        if supplied_as := clean_value(item.get("supplied_as")):
            add_amount_quantity(output, supplied_as)

        if passage := clean_value(item.get("passage_history")):
            passage_text = passage if "passage" in passage.casefold() else f"{passage} passage"
            insert_value(output, "sampleState", passage_text, extend=True)

        if mutations := clean_value(item.get("mutations")):
            insert_value(output, "associatedGenotype", mutations)

        if vector := clean_value(item.get("vector")):
            insert_value(output, "isBasedOn", {"@type": "BioChemEntity", "name": vector})

        if immunogen := clean_value(item.get("immunogen")):
            insert_value(output, "isBasedOn", {"@type": "BioChemEntity", "name": immunogen})

        if produced_in := clean_value(item.get("produced_in")):
            insert_value(output, "species", {"name": produced_in})

        if compatibility := clean_value(item.get("compatibility")):
            insert_value(output, "experimentalPurpose", compatibility)

        if comments := clean_value(item.get("comments")):
            insert_value(output, "description", comments, extend=True)

        if concentration := clean_value(item.get("concentration")):
            add_concentration_quantity(output, concentration)

        quantity_available = maybe_number(item.get("quantity_available"))
        if quantity_available is not None:
            add_sample_quantity(output, "# Aliquots Available", value=quantity_available, unit_text="Aliquots")
            if quantity_available >= 1:
                output["creativeWorkStatus"] = "Available"
            elif quantity_available == 0:
                output["creativeWorkStatus"] = "Backordered"

        if quantity_minimum := maybe_number(item.get("quantity_minimum")):
            add_sample_quantity(output, "minimum order", value=quantity_minimum, unit_text="minimum order")

        if purification := expand_abbreviation(item, "purification", item.get("purification")):
            insert_value(output, "sampleProcess", f"{purification} purification", extend=True)

        if antibody_type := expand_abbreviation(item, "antibody_type", item.get("antibody_type")):
            insert_value(output, "keywords", antibody_type)

        if specificity := clean_value(item.get("specificity")):
            insert_value(output, "keywords", specificity)

        if protein := clean_value(item.get("protein")):
            insert_value(output, "hasPart", {"@type": "BioChemEntity", "name": protein})

        if segment := clean_value(item.get("segment")):
            insert_value(output, "hasPart", {"@type": "BioChemEntity", "name": segment})

    def add_publications(self, output, item):
        for pmid in split_pmids(item.get("publication_pmid")):
            insert_value(
                output,
                "isBasedOn",
                {
                    "@type": "ScholarlyArticle",
                    "identifier": f"PMID:{pmid}",
                    "pmid": pmid,
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                },
            )

    def add_bei_link(self, output, bei_number):
        if not bei_number:
            return
        bei_url = f"https://www.beiresources.org/Catalog/BEIResources/{bei_number}.aspx"
        output["sameAs"] = bei_url
        insert_value(
            output,
            "isRelatedTo",
            {
                "@type": "Sample",
                "identifier": self.bei_overlap_id(bei_number),
                "name": f"BEI {bei_number}",
                "url": bei_url,
                "includedInDataCatalog": {
                    "@type": "DataCatalog",
                    "name": "BEI Resources",
                    "url": self.BEI_HOME_URL,
                },
                "relationship": "also available from",
            },
        )

    def add_gene_segments(self, output, item):
        for field in SEGMENT_FIELDS:
            value = clean_value(item.get(field))
            if not value:
                continue
            insert_value(
                output,
                "isBasedOn",
                {"@type": "BioChemEntity", "name": value, "relationship": f"{field} segment source"},
            )
