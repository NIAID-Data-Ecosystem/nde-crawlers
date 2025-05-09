import datetime
import logging

from dateutil.parser import parse as date_parse

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(name)s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("nde-logger")


def parse(record, dataset_name, _id, url):
    output = {
        "@context": "http://schema.org/",
        "@type": "Dataset",
        "url": url,
        "_id": _id,
        "identifier": _id,
        "includedInDataCatalog": {
            "@type": "DataCatalog",
            "name": "Omics Discovery Index (OmicsDI)",
            "url": "https://www.omicsdi.org/",
            "versionDate": datetime.date.today().isoformat(),
            "dataset": url,
        },
        "distribution": {
            "@type": "DataDownload",
            "encodingFormat": "XML",
            "contentUrl": f"https://www.omicsdi.org/ws/dataset/{dataset_name}/{_id}",
        },
    }

    keywords_set = set()
    authors = []
    mt = []
    vm = []
    species = []
    ia = []
    hc = []
    funding_list = []
    citations = []
    spatial_coverage = []
    is_related_to = []

    sd_publisher = {}
    if name := record.get("database"):
        sd_publisher["name"] = name
        sd_publisher["url"] = f"https://www.omicsdi.org/dataset/{dataset_name}/{_id}"
        sd_publisher["@type"] = "Organization"
        output["sdPublisher"] = sd_publisher

    if description := record.get("description"):
        output["description"] = description

    if name := record.get("name"):
        output["name"] = name

    if dates := record.get("dates"):
        date_fields = {
            "publication": "datePublished",
            "submission": "dateCreated",
            "created": "dateCreated",
            "creation": "dateCreated",
            "last_modification": "dateModified",
            "last_modified": "dateModified",
            "modification": "dateModified",
            "updated": "dateModified",
            "release": "datePublished",
        }

        for key, output_field in date_fields.items():
            if date_value := dates.get(key):
                try:
                    parsed_date = (
                        date_parse(date_value, ignoretz=True).date().isoformat()
                    )
                    output[output_field] = parsed_date
                except Exception as e:
                    logger.error(f"Error parsing date for {key}: {e}")

    if additional := record.get("additional"):

        if isprivate := additional.get("isPrivate"):
            if isprivate == "true":
                return None

        # List of keys to check in `additional`
        keyword_keys = [
            "submitter_keywords",
            "curator_keywords",
            "instrument_platform",
            "modification",
            "software",
            "keywords",
            "metabolite_name",
            "PerturbationType_Level1",
            "PerturbationType_Level2",
            "PerturbationType_Level3",
            "keyword",
        ]

        # Iterate over the keys and update `keywords_set`
        for key in keyword_keys:
            if keywords := additional.get(key):
                keywords_set.update(keywords)

        def parse_authors(author_names, affiliation=None, orcid=None):
            authors = []
            for i, author_name in enumerate(author_names):
                author = {"@type": "Person"}
                author["name"] = author_name
                if affiliation:
                    author["affiliation"] = {
                        "name": (
                            affiliation[i]
                            if len(affiliation) == len(author_names)
                            else affiliation[0]
                        )
                    }
                if orcid:
                    author["identifier"] = (
                        f"https://orcid.org/{orcid[i]}"
                        if len(orcid) == len(author_names)
                        else f"https://orcid.org/{orcid[0]}"
                    )
                authors.append(author)
            return authors

        if author_names := additional.get("submitter"):
            orcid = None
            if record.get("cross_references") and record.get("cross_references").get(
                "ORCID"
            ):
                orcid = record.get("cross_references").get("ORCID")

            affiliation = additional.get("submitter_affiliation")
            authors.extend(parse_authors(author_names, affiliation, orcid))

        if author_names := additional.get("labhead"):
            affiliation = additional.get("labhead_affiliation")
            authors.extend(parse_authors(author_names, affiliation))

        if author_names := additional.get("contact_person"):
            authors.extend(parse_authors(author_names))

        if author_names := additional.get("corresponding_author_name"):
            authors.extend(parse_authors(author_names))

        if author_names := additional.get("first_author"):
            authors.extend(parse_authors(author_names))

        if not authors:
            if author_names := additional.get("author"):
                authors.extend(parse_authors(author_names))

        # List of keys to check in `additional`
        mt_keys = [
            "quantification_method",
            "technology_type",
            "dataset_type",
            "study_design",
            "study_type",
        ]

        # Iterate over the keys and append to `mt`
        for key in mt_keys:
            if mt_names := additional.get(key):
                for mt_name in mt_names:
                    mt.append({"name": mt_name})

        protocols = {
            "capillary_electrophoresis_protocol": {
                "url": "http://purl.obolibrary.org/obo/CHMO_0000702",
                "name": "capillary electrophoresis-mass spectrometry",
            },
            "chromatography_protocol": {
                "url": "http://purl.obolibrary.org/obo/MMO_0000288",
                "name": "chromatography",
            },
            "hitology_protocol": {
                "url": "http://purl.obolibrary.org/obo/MMO_0000496",
                "name": "histology",
            },
            "in_vivo_magnetic_resonance_assay_protocol": {
                "url": "http://purl.obolibrary.org/obo/MMO_0000019",
                "name": "magnetic resonance imaging",
            },
            "in_vivo_magnetic_resonance_spectroscopy_protocol": {
                "url": "http://purl.obolibrary.org/obo/NCIT_C16810",
                "name": "Magnetic Resonance Spectroscopy",
            },
            "lipidomics_uplc_ms_protocol": {
                "url": "http://purl.obolibrary.org/obo/MMO_0000721",
                "name": "ultra-performance liquid chromatography-mass spectrometry",
            },
            "magnetic_resonance_imaging_protocol": {
                "url": "http://purl.obolibrary.org/obo/MMO_0000019",
                "name": "magnetic resonance imaging",
            },
            "mass_spectronomy_protocol": {
                "url": "http://purl.obolibrary.org/obo/MMO_0000534",
                "name": "mass spectrometry",
            },
            "mass_spec_protocol": {
                "url": "http://purl.obolibrary.org/obo/MMO_0000534",
                "name": "mass spectrometry",
            },
            "matrix_assisted_laser_desorption_ionization_protocol": {
                "url": "http://purl.obolibrary.org/obo/CHMO_0000054",
                "name": "matrix-assisted laser desorption-ionisation mass spectrometry",
            },
            "nmr_assay_protocol": {
                "url": "http://purl.obolibrary.org/obo/MMO_0000709",
                "name": "nuclear magnetic resonance spectroscopy",
            },
            "nmr_spec_protocol": {
                "url": "http://purl.obolibrary.org/obo/MMO_0000709",
                "name": "nuclear magnetic resonance spectroscopy",
            },
        }

        # Iterate over the lookup dictionary
        for key, value in protocols.items():
            if additional.get(key):  # Check if the key exists in `additional`
                mt.append(value)

        # Define a lookup dictionary for prefixes
        imaging_protocols = {
            "3d_desi_imaging": {
                "url": "http://purl.obolibrary.org/obo/CHMO_0000484",
                "name": "desorption electrospray ionisation mass spectrometry",
            },
            "3d_maldi_imaging": {
                "url": "http://purl.obolibrary.org/obo/CHMO_0000054",
                "name": "matrix-assisted laser desorption-ionisation imaging mass spectrometry",
            },
        }

        # Iterate over the lookup dictionary
        for prefix, value in imaging_protocols.items():
            if any(
                key.startswith(prefix) and additional.get(key)
                for key in additional.keys()
            ):
                mt.append(value)
                break

        if sameAs := additional.get("full_dataset_link"):
            output["sameAs"] = sameAs

        if vm_names := additional.get("omics_type"):
            for vm_name in vm_names:
                vm.append({"name": vm_name})

        if vm_names := additional.get("study_factor"):
            for vm_name in vm_names:
                vm.append({"name": vm_name})

        for key in ["species", "Organism", "taxonomy"]:
            if species_names := additional.get(key):
                for species_name in species_names:
                    species.append({"name": species_name})

        if additional.get("GenotypeID") and not additional.get("Organism"):
            for infectious_agent in additional.get("GenotypeID"):
                ia.append({"name": infectious_agent})
        elif species_names := additional.get("GenotypeID"):
            for species_name in species_names:
                species.append({"name": species_name})

        if (
            additional.get("GenotypeID_Systematic")
            and not additional.get("Organism")
            and not additional.get("GenotypeID")
        ):
            for infectious_agent in additional.get("GenotypeID_Systematic"):
                ia.append({"name": infectious_agent})
        elif species_names := additional.get("GenotypeID_Systematic"):
            for species_name in species_names:
                species.append({"name": species_name})

        for key in ["condition", "disease"]:
            if hc_names := additional.get(key):
                for hc_name in hc_names:
                    hc.append({"name": hc_name})

        if description := additional.get("description"):
            if output.get("description"):
                output["description"] += description[0]
            else:
                output["description"] = description[0]

        if legend := additional.get("legend"):
            if not output.get("description"):
                output["description"] = legend[0]

        if doi := additional.get("doi"):
            output["doi"] = doi

        if funder_names := additional.get("funding"):
            if identifiers := additional.get("funding_grand_id"):
                if len(identifiers) == len(funder_names):
                    pass
                else:
                    identifiers = None
            else:
                identifiers = None

            for i, funder_name in enumerate(funder_names):
                funding = {"funder": {"name": funder_name}}
                if identifiers:
                    funding["identifier"] = identifiers[i]
                funding_list.append(funding)

        if locations := additional.get("location"):
            for location in locations:
                spatial_coverage.append({"name": location})

        if dois := additional.get("publication_doi"):
            for doi in dois:
                citation = {"doi": doi}
                citations.append(citation)
        if pmids := additional.get("publication_pubmed"):
            for pmid in pmids:
                citation = {"pmid": pmid}
                citations.append(citation)

        if pmids := additional.get("pubmed"):
            for pmid in pmids:
                citation = {"pmid": pmid}
                citations.append(citation)

        if start_date := additional.get("study_start_year"):
            try:
                start_date = date_parse(start_date[0], ignoretz=True).date().isoformat()
                output["temporalCoverage"] = {
                    "temporalInterval": {"startDate": start_date}
                }
            except Exception as e:
                logger.error(f"Error parsing date for study_start_year: {e}")

    if cross_references := record.get("cross_references"):
        if pmids := cross_references.get("pubmed"):
            output["pmids"] = ",".join(pmids)

        if doi := cross_references.get("doi"):
            is_related_to.append({"doi": doi})

        is_related_to_keys = [
            "Chinese Clinical Trial Register",
            "Chinese Medicine Clinical Trials Registry",
            "cl",
            "Clinical Research Information Service",
            "Clinical Trials Registry - India",
            "ClinicalTrials___gov",
            "Deutschen Register Klinischer Studien",
            "EU Clinical Trials Register",
            "International Clinical Trials Registry Platform",
            "Iranian Registry of Clinical Trials",
            "Primary Registries Network",
            "Japan Registry of Clinical Trials",
            "nct",
            "Netherlands National Trial Register",
            "Registro Brasileiro de Ensaios Clínicos",
            "Thai Clinical Trials Register",
        ]

        # Iterate over the keys and append identifiers
        for key in is_related_to_keys:
            if identifiers := cross_references.get(key):
                for identifier in identifiers:
                    is_related_to.append({"identifier": identifier})

        # Define a lookup dictionary for cross_references
        is_related_to_lookup = {
            "arrayexpress": "https://www.omicsdi.org/dataset/biostudies-arrayexpress/{}",
            "biomodels__db": "https://www.omicsdi.org/dataset/biomodels/{}",
            "bioproject": "https://www.omicsdi.org/dataset/project/{}",
            "massive": "https://www.omicsdi.org/dataset/massive/{}",
            "Massive": "https://www.omicsdi.org/dataset/massive/{}",
            "Pride": "https://www.omicsdi.org/dataset/pride/{}",
            "PrideTest": "https://www.omicsdi.org/dataset/pride/{}",
            "PrideArchive": "https://www.omicsdi.org/dataset/pride/{}",
            "PeptideAtlas": "https://www.omicsdi.org/dataset/peptide_atlas/{}",
            "Metabolights": "https://www.omicsdi.org/dataset/metabolights_dataset/{}",
            "EGA": "https://www.omicsdi.org/dataset/ega/{}",
            "NODE": "https://www.omicsdi.org/dataset/node/{}",
        }

        # Iterate over the lookup dictionary
        for key, url_template in is_related_to_lookup.items():
            if identifiers := cross_references.get(key):
                for identifier in identifiers:
                    is_related_to.append(
                        {
                            "identifier": identifier,
                            "url": url_template.format(identifier),
                        }
                    )

        if hc_names := cross_references.get("hp"):
            for hc_name in hc_names:
                hc.append({"identifier": hc_name})

    # add lists to output
    if keywords_set:
        keywords_list = list(keywords_set)
        if len(keywords_list) > 1:
            output["keywords"] = keywords_list
        else:
            output["keywords"] = keywords_list[0]

    # Map variables to their corresponding output keys
    output_mappings = {
        "author": authors,
        "measurementTechnique": mt,
        "variableMeasured": vm,
        "species": species,
        "infectiousAgent": ia,
        "healthCondition": hc,
        "funding": funding_list,
        "citation": citations,
        "spatialCoverage": spatial_coverage,
        "isRelatedTo": is_related_to,
    }

    # Iterate over the mappings and add non-empty values to the output
    for output_key, value in output_mappings.items():
        if value:  # Check if the value is not empty
            output[output_key] = value

    return output
