import datetime
import json
import logging
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from typing import Any, Dict, Generator, List, Optional

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("bioproject-parser")

BIOPROJECT_XML_URL = "https://ftp.ncbi.nlm.nih.gov/bioproject/bioproject.xml"
BIOPROJECT_XML_PATH = "bioproject.xml"
BIOPROJECT_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/bioproject/summary.txt"

# Global variable to cache summary data
_summary_data_cache = None

def load_summary_data() -> Dict[str, Dict[str, str]]:
    """Load summary data as a lookup table for enrichment."""
    global _summary_data_cache

    if _summary_data_cache is not None:
        return _summary_data_cache

    logger.info("Loading BioProject summary data from %s", BIOPROJECT_SUMMARY_URL)
    summary_data = {}

    try:
        with urllib.request.urlopen(BIOPROJECT_SUMMARY_URL) as response:
            lines = response.read().decode('utf-8').splitlines()

        # Skip header line
        for line in lines[1:]:
            parts = line.split('\t')
            if len(parts) >= 7:
                accession = parts[2]  # Project Accession
                summary_data[accession] = {
                    'organism_name': parts[0],
                    'taxid': parts[1],
                    'project_id': parts[3],
                    'project_type': parts[4],
                    'data_type': parts[5],
                    'date': parts[6]
                }

        logger.info("Loaded summary data for %d BioProject records", len(summary_data))
        _summary_data_cache = summary_data
        return summary_data

    except Exception as e:
        logger.warning("Failed to load summary data: %s", e)
        return {}

def enrich_with_summary_data(output: Dict[str, Any], accession: str, summary_data: Dict[str, Dict[str, str]]) -> None:
    """Enrich output with summary data if available and fields are missing."""
    if not accession or accession not in summary_data:
        return

    summary = summary_data[accession]
    enriched_fields = []

    # Fill in missing organism/species data
    if 'species' not in output and summary.get('organism_name') and summary.get('taxid'):
        output['species'] = {
            'name': summary['organism_name'],
            'identifier': f"NCBITaxon:{summary['taxid']}"
        }
        enriched_fields.append('species')

    # Fill in missing measurement technique from project type
    if 'measurementTechnique' not in output and summary.get('project_type'):
        project_type = summary['project_type'].lstrip('e').strip()
        if project_type and project_type != 'Primary submission':
            output['measurementTechnique'] = {'name': project_type}
            enriched_fields.append('measurementTechnique')

    # Fill in missing variable measured from data type
    if 'variableMeasured' not in output and summary.get('data_type'):
        data_type = summary['data_type'].strip()
        if data_type:
            output['variableMeasured'] = {'name': data_type}
            enriched_fields.append('variableMeasured')

    # Fill in missing creation date
    if 'datePublished' not in output and summary.get('date'):
        try:
            # Convert YYYY/MM/DD to YYYY-MM-DD
            date_str = summary['date'].replace('/', '-')
            # Validate the date format
            datetime.datetime.strptime(date_str, '%Y-%m-%d')
            output['datePublished'] = date_str
            enriched_fields.append('datePublished')
        except ValueError:
            pass  # Skip invalid dates

    if enriched_fields:
        logger.debug("Enriched %s with fields: %s", accession, ', '.join(enriched_fields))


def download_bioproject_xml():
    """Download the BioProject XML file from NCBI FTP."""
    logger.info("Downloading BioProject XML from %s", BIOPROJECT_XML_URL)
    urllib.request.urlretrieve(BIOPROJECT_XML_URL, BIOPROJECT_XML_PATH)
    logger.info("Downloaded BioProject XML file to %s", BIOPROJECT_XML_PATH)


def get_ids():
    """Get iterator for parsing BioProject XML."""
    logger.info("Streaming records from BioProject XML")
    return ET.iterparse(BIOPROJECT_XML_PATH, events=("end",))


def get_related_sra_projects(accession: str) -> List[str]:
    """
    Query NCBI Entrez API to find related SRA projects for a BioProject accession.
    """
    try:
        # First, search for SRA studies linked to this BioProject
        search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        search_params = {
            'db': 'sra',
            'term': f'{accession}[BioProject]',
            'retmax': '1000',  # Adjust as needed
            'retmode': 'xml'
        }

        search_query = search_url + '?' + urllib.parse.urlencode(search_params)

        with urllib.request.urlopen(search_query) as response:
            search_result = response.read()

        search_root = ET.fromstring(search_result)

        # Extract the SRA UIDs
        uids = [uid.text for uid in search_root.findall('.//Id')]

        if not uids:
            logger.info("No SRA projects found for BioProject %s", accession)
            return []

        # Now get the SRA study accessions using esummary
        summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        summary_params = {
            'db': 'sra',
            'id': ','.join(uids),
            'retmode': 'xml'
        }

        summary_query = summary_url + '?' + urllib.parse.urlencode(summary_params)

        with urllib.request.urlopen(summary_query) as response:
            summary_result = response.read()

        summary_root = ET.fromstring(summary_result)

        # Extract SRP accessions from the summary
        srp_accessions = []
        for doc_sum in summary_root.findall('.//DocSum'):
            for item in doc_sum.findall('.//Item'):
                if item.get('Name') == 'ExpXml':
                    # Parse the ExpXml to find Study accession
                    exp_xml = item.text
                    if exp_xml:
                        try:
                            exp_root = ET.fromstring(exp_xml)
                            study_acc = exp_root.findtext('.//Study/@acc')
                            if study_acc and study_acc.startswith('SRP'):
                                srp_accessions.append(study_acc)
                        except ET.ParseError:
                            continue

        # Remove duplicates
        logger.info("Found %d related SRA projects for BioProject %s", len(srp_accessions), accession)
        return list(set(srp_accessions))

    except Exception as e:
        logger.warning(f"Failed to get related SRA projects for {accession}: {e}")
        return []


def debug_xml_structure():
    """Debug function to inspect the XML structure."""
    try:
        # Parse just the beginning of the file to see structure
        with open(BIOPROJECT_XML_PATH, 'rb') as f:
            # Read first few KB to see the structure
            sample = f.read(5000).decode('utf-8', errors='ignore')
            print("First 5000 characters of XML:")
            print(sample)
            print("\n" + "="*50 + "\n")

        # Try to parse and see what the root element actually is
        tree = ET.parse(BIOPROJECT_XML_PATH)
        root = tree.getroot()
        print(f"Root element tag: {root.tag}")
        print(f"Root element attributes: {root.attrib}")

        # Print first few child elements
        print("First 5 child elements:")
        for i, child in enumerate(root):
            if i >= 5:
                break
            print(f"  {i}: {child.tag} - {child.attrib}")

    except Exception as e:
        print(f"Error debugging XML: {e}")


def parse() -> Generator[Dict[str, Any], None, None]:
    """Parse BioProject XML and yield structured records."""
    # Download the XML file
    download_bioproject_xml()

    # Load summary data for enrichment
    summary_data = load_summary_data()

    context = get_ids()
    logger.info("Starting to parse BioProject XML")
    _, root = next(context)
    count = 0

    for event, elem in context:
        if elem.tag != "Package":
            continue

        # Find the nested Project element within the Package
        project_elem = elem.find(".//Project/Project")
        if project_elem is None:
            continue

        count += 1
        if count % 50 == 0:
            logger.info("Parsed %s records", count)

        output = {
            "@type": "Dataset",
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "NCBI BioProject",
                "url": "https://www.ncbi.nlm.nih.gov/bioproject/",
                "versionDate": datetime.date.today().isoformat(),
            }
        }

        # Extract accession and build identifiers
        archive_id = project_elem.find(".//ProjectID/ArchiveID")
        accession = archive_id.get("accession") if archive_id is not None else None

        if accession:
            output["identifier"] = accession
            output["_id"] = f"ncbi_bioproject_{accession}"
            output["url"] = f"https://www.ncbi.nlm.nih.gov/bioproject/{accession}"
            output["includedInDataCatalog"]["dataset"] = output["url"]

        # Extract title and name
        title = project_elem.findtext(".//ProjectDescr/Title")
        name = project_elem.findtext(".//ProjectDescr/Name")
        if title and name:
            output["name"] = f"{title}. {name}"
        else:
            output["name"] = title or name

        # Extract description
        if desc := project_elem.findtext(".//ProjectDescr/Description"):
            # Clean up HTML tags and formatting
            cleaned_desc = desc.strip()
            # Remove HTML tags
            cleaned_desc = re.sub(r'<[^>]+>', '', cleaned_desc)
            # Replace multiple whitespace with single space
            cleaned_desc = re.sub(r'\s+', ' ', cleaned_desc)
            # Remove carriage returns and line feeds
            cleaned_desc = cleaned_desc.replace('\r', '').replace('\n', ' ')
            # Final cleanup of extra spaces
            cleaned_desc = cleaned_desc.strip()

            if cleaned_desc:
                output["description"] = cleaned_desc

        # Extract dates
        submission_elem = project_elem.find(".//Submission")
        if submission_elem is not None:
            if submitted_date := submission_elem.get("submitted"):
                output["datePublished"] = submitted_date
            if last_update_date := submission_elem.get("last_update"):
                output["dateModified"] = last_update_date

        # Extract data types (variableMeasured)
        # Priority: ProjectDataTypeSet.DataType > IntendedDataTypeSet.DataType > Objectives.Data.@data_type
        data_types = []

        # First try ProjectDataTypeSet.DataType
        data_types = [dt.text.strip() for dt in project_elem.findall(".//ProjectDataTypeSet/DataType") if dt.text]

        # Fallback to IntendedDataTypeSet.DataType
        if not data_types:
            data_types = [dt.text.strip() for dt in project_elem.findall(".//IntendedDataTypeSet/DataType") if dt.text]

        # Final fallback to Objectives.Data.@data_type
        if not data_types:
            data_types = [
                dt.get("data_type").lstrip("e").strip()
                for dt in project_elem.findall(".//Objectives/Data")
                if dt.get("data_type")
            ]

        if data_types:
            # Convert to objects with name property
            variable_measured = [{"name": dt} for dt in data_types]
            output["variableMeasured"] = variable_measured if len(variable_measured) > 1 else variable_measured[0]

        # Extract measurement technique
        method = project_elem.find(".//ProjectTypeSubmission/Method")
        if method is not None:
            if method_type := method.get("method_type"):
                output["measurementTechnique"] = {"name": method_type.lstrip("e")}

        # Extract species/organism information
        organism = project_elem.find(".//ProjectTypeSubmission/Target/Organism")
        if organism is not None:
            species = {}
            if taxid := organism.get("taxID"):
                species["identifier"] = f"NCBITaxon:{taxid}"
            elif species_id := organism.get("species"):
                species["identifier"] = f"NCBITaxon:{species_id}"

            if org_name := organism.findtext("OrganismName"):
                species["name"] = org_name

            if species:
                output["species"] = species

        # Extract funding information
        grant = project_elem.find(".//ProjectDescr/Grant")
        if grant is not None:
            funding = {}
            if grant_id := grant.get("GrantId"):
                funding["identifier"] = grant_id
                logger.debug(f"Found grant ID: {grant_id}")
            if gtitle := grant.findtext("Title"):
                funding["name"] = gtitle
                logger.debug(f"Found grant title: {gtitle}")

            agency = grant.find("Agency")
            if agency is not None:
                funder = {}
                if abbr := agency.get("abbr"):
                    funder["alternateName"] = abbr
                    logger.debug(f"Found agency abbr: {abbr}")
                if agency.text:
                    funder["name"] = agency.text
                    logger.debug(f"Found agency name: {agency.text}")
                if funder:
                    funding["funder"] = funder

            if funding:
                output["funding"] = funding
                logger.debug(f"Added funding to output: {funding}")
            else:
                logger.debug("No funding data found - empty funding dict")
        else:
            logger.debug("No Grant element found")

        # Extract authors
        author_list = []

        # Add PI from grant
        if grant is not None:
            pi = grant.find("PI")
            if pi is not None:
                pi_author = {"@type": "Person"}
                if given := pi.findtext("Given"):
                    pi_author["givenName"] = given
                if last := pi.findtext("Last"):
                    pi_author["familyName"] = last
                if affil := pi.get("affil"):
                    pi_author["affiliation"] = {"name": affil}
                if len(pi_author) > 1:  # More than just @type
                    author_list.append(pi_author)

        # Also check Provider field for PI information
        target = project_elem.find(".//ProjectTypeSubmission/Target")
        if target is not None:
            provider = target.findtext("Provider")
            if provider and provider.startswith("PI:"):
                pi_author = {
                    "@type": "Person",
                    "name": provider.replace("PI:", "").strip()
                }
                author_list.append(pi_author)

        # Add organization
        org = project_elem.find(".//Submission/Description/Organization")
        if org is not None:
            org_author = {"@type": "Organization"}

            # Get organization name
            org_name = org.findtext("Name")
            if not org_name:
                # Try Name/@abbr
                name_elem = org.find("Name")
                if name_elem is not None:
                    org_name = name_elem.get("abbr")
            if not org_name:
                # Try Name/#text
                name_elem = org.find("Name")
                if name_elem is not None and name_elem.text:
                    org_name = name_elem.text

            if org_name:
                org_author["name"] = org_name

                if role := org.get("role"):
                    org_author["role"] = role

                author_list.append(org_author)

        if author_list:
            output["author"] = author_list

        # Extract publications
        publications = []
        for pub in project_elem.findall(".//ProjectDescr/Publication"):
            pub_entry = {}
            if pub_id := pub.get("id"):
                if db := pub.get("DbType"):
                    if db == "ePubmed":
                        pub_entry["pmid"] = pub_id
                    elif db == "ePubMedCentral":
                        pub_entry["pmcid"] = pub_id

            if pub_entry:
                publications.append(pub_entry)

        if publications:
            output["citation"] = publications

        # Extract hierarchical relationships and related SRA projects
        numeric_id = project_elem.findtext(".//ProjectID/ArchiveID[@id]")
        has_parts = []
        part_of = []

        # Process BioProject hierarchical relationships
        for link in project_elem.findall(".//ProjectLinks/Link"):
            ref_acc = link.findtext("ProjectIDRef[@accession]")
            hierarchical = link.find("Hierarchical")

            if hierarchical is not None and ref_acc:
                if hierarchical.get("type") == "TopAdmin":
                    member = hierarchical.find("MemberID")
                    if member is not None:
                        member_id = member.get("id")
                        if member_id == numeric_id:
                            # This project is the parent of the referenced project
                            has_parts.append({"identifier": ref_acc})
                        else:
                            # This project is a child of the referenced project
                            part_of.append({"identifier": ref_acc})

        # Add related SRA projects to hasPart (Option 3 implementation)
        # if accession:
        #     related_sra_projects = get_related_sra_projects(accession)
        #     if related_sra_projects:
        #         # Add SRA project IDs to hasPart
        #         # Format them as portal URLs or identifiers as needed
        #         sra_parts = [{"identifier": f"ncbi_sra_{srp_id}"} for srp_id in related_sra_projects]
        #         has_parts.extend(sra_parts)

        if has_parts:
            output["hasPart"] = has_parts
        if part_of:
            output["isPartOf"] = part_of

        # Extract access conditions
        access = project_elem.findtext(".//Submission/Description/Access")
        if access:
            if access.lower() == "public":
                output["conditionsOfAccess"] = "Open"

        # Extract keywords from relevance and properties
        keywords = []

        # Medical relevance
        if project_elem.findtext(".//ProjectDescr/Relevance/Medical") == "yes":
            keywords.append("medical")

        # Sample scope as keyword
        target_elem = project_elem.find(".//ProjectTypeSubmission/Target")
        if target_elem is not None:
            sample_scope = target_elem.get("sample_scope")
            if sample_scope and sample_scope.strip():
                keywords.append(sample_scope.lstrip("e").strip().lower())

        if keywords:
            output["keywords"] = keywords

        # Enrich with summary data for any missing fields
        if accession and summary_data:
            enrich_with_summary_data(output, accession, summary_data)

        yield output
        elem.clear()
        root.clear()

    logger.info("Finished parsing %s BioProject records", count)
