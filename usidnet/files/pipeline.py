# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: https://docs.scrapy.org/en/latest/topics/item-pipeline.html


# useful for handling different item types with a single interface
# mapping https://docs.google.com/spreadsheets/d/19hQf4sQ6ZLwvt8ADYyvooYrqmQCsxh386kYnycX3xmg/edit#gid=0

import datetime
import logging
import re

import dateutil.parser

logger = logging.getLogger("nde-logger")


__all__ = [
    "USIDNETItemProcessorPipeline",
]


def insert_value(d, key, value, extend=False):
    """Insert a value into a dictionary, handling existing keys by converting to lists or extending strings as needed."""

    if key in d and not extend:
        if isinstance(d[key], list):
            if isinstance(value, list):
                for item in value:
                    if item not in d[key]:
                        d[key].append(item)
            elif value not in d[key]:
                d[key].append(value)
        else:
            if isinstance(value, list):
                d[key] = [d[key]] + [v for v in value if v != d[key]]
            elif d[key] != value:
                d[key] = [d[key], value]
    elif d.get(key) and extend:
        d[key] = (d.get(key) + " " + value).strip()
    else:
        d[key] = value


def _to_iso_date(val):
    if val is None:
        return None
    try:
        dt = dateutil.parser.parse(val, ignoretz=True).date().isoformat()
    except (dateutil.parser.ParserError, TypeError):
        logger.warning(f"Could not parse date: {val}")
        return None
    return dt


# Values that mean "no condition / not applicable / unknown" — used to gate
# inclusion of healthCondition / associatedGenotype values per the mapping notes.
_EXCLUSION_TOKENS = frozenset(
    {
        "",
        "absent",
        "no",
        "no data",
        "none",
        "none known",
        "no notes:",
        "not applicable",
        "not done",
        "not examined",
        "not given",
        "not tested",
        "not tried/not given",
        "not yet",
        "not yet reported",
        "n/a",
        "na",
        "null",
        "unknown",
        "unspecified",
    }
)

def _normalize(value):
    if value is None:
        return ""
    return str(value).strip()


def _clean_token(value):
    v = _normalize(value).casefold()
    # tolerate trailing "notes: ..." text scrapy sometimes captures
    v = v.split("notes:", 1)[0].strip()
    return v


def _is_excluded(value):
    return _clean_token(value) in _EXCLUSION_TOKENS


# healthCondition columns where the COLUMN NAME is the disease and the VALUE is
# a presence indicator (e.g. "Epilepsy": "present"). Default behavior for
# anything mapped to healthCondition.name unless listed in _HC_VALUE_IS_DISEASE.
_HC_VALUE_IS_DISEASE = frozenset(
    {
        "Cerebrovascular diagnosis",
        "diagnosis",
        "Ischemic stroke subtype",
        "ischemic stroke subtype based on TOAST",
        "Known Genetic Syndrome",
        "Known genetic syndrome",
        "Neurological Symptoms",
        "Non-dystonia syndrome",
        "Other dystonia syndrome",
        "Parkinsonism clinical diagnosis",
        "Primary clinical diagnosis",
        "Primary Dystonia type",
        "Secondary neurological diagnoses",
    }
)

_HEALTH_CONDITION_FIELDS = frozenset(
    {
        "ALS/other MND",
        "Activation tremor",
        "Age-related macular degeneration (AMD)",
        "Alzheimer disease",
        "Alzheimer's disease",
        "Amyotrophic lateral sclerosis",
        "Anemia",
        "Anxiety",
        "Arthritis",
        "Arthritis, unspecified",
        "Asthma",
        "Ataxia",
        "Atherosclerosis",
        "Atonic",
        "Atrial fibrillation",
        "Atrial septal defect",
        "Attention Deficit Hyperactivity Disorder (ADHD)",
        "Attention Deficit/Hyperactivity Disorder",
        "Atypical Dysphonia",
        "Autism",
        "Benign prostate hypertrophy",
        "Bipolar (manic-depressive)",
        "Bradykinesis",
        "Brain aneurysm",
        "Brain tumor",
        "Cancer",
        "Cataracts",
        "Cerebellar signs (other than activation tremor)",
        "Cerebrovascular Disorders, including stroke and aneurysm",
        "Cerebrovascular diagnosis",
        "Chronic Multiple Tics",
        "Chronic obstructive pulmonary disease (COPD)",
        "Clonic",
        "Complex febrile seizure",
        "Congestive heart disease",
        "Corticobasal ganglionic degeneration",
        "Cysts",
        "Degenerative joint disease",
        "Delusions",
        "Dementia",
        "Depression",
        "Diabetes",
        "Diabetes mellitus",
        "Diabetic retinopathy",
        "Diffuse Lewy Body disease",
        "Dystonia",
        "Early severe dementia",
        "Emotional incontinence",
        "Emphysema",
        "Epilepsy",
        "Fluctuations in attention or alertness",
        "Focal Dystonia",
        "Focal dyscognitive",
        "Focal evolving to generalized convulsion",
        "Focal motor /autonomic/aura",
        "Focal sensory deficit",
        "Foot Dystonia present",
        "Foot tremor",
        "Foot tremor type",
        "Frailty",
        "Frontal lobe dysfunction (bradyphrenia, forced grasping, perseveration, utilization behavior)",
        "Gait difficulties",
        "Gastritis",
        "Gastroesophageal reflux disease",
        "Gaze palsy",
        "Glaucoma",
        "Hand Dystonia present",
        "Hand tremor present",
        "Hand tremor type",
        "Hearing loss",
        "Heart attack",
        "Heart disease",
        "Hemi Dystonia",
        "High blood pressure",
        "High cholesterol",
        "Hydrocephalus",
        "Hypercholesterolemia",
        "Hypertension",
        "Hysterectomy",
        "If yes, then type:",
        "Infantile spasms",
        "Insomnia",
        "Ischemic stroke subtype",
        "Jaw Dystonia present",
        "Jaw tremor present",
        "Jaw tremor type",
        "Kidney/bladder stones",
        "Known Genetic Syndrome",
        "Known genetic syndrome",
        "Laryngeal tremor",
        "Larynx Dystonia present",
        "Lower Face Dystonia present",
        "Lower Face Tremor",
        "Lower Face tremor type",
        "Macular degeneration",
        "Memory loss",
        "Migraine",
        "Migraines",
        "Motor Neuron Disorders, including ALS",
        "Multiple myeloma",
        "Multiple sclerosis",
        "Multiple system atrophy",
        "Muscle disease",
        "Muscular tension Dysphonia",
        "Myocardial infarction",
        "Myoclonic",
        "Myoclonus",
        "Myopia",
        "Neck Dystonia present",
        "Neck tremor present",
        "Neck tremor type",
        "Neuroleptic treatment onset",
        "Neurological Symptoms",
        "Non-dystonia syndrome",
        "Nonepileptiform abnormalities focal slowing",
        "Obesity",
        "Obsessive Compulsive",
        "Obsessive Compulsive Disorder (OCD)",
        "Obsessive compulsive disorder",
        "Oculogyric crisis",
        "Osteoarthritis",
        "Osteoporosis",
        "Other Tic disorder",
        "Other dementia",
        "Other dystonia syndrome",
        "Other neurodegenerative disease",
        "Other types of Psychiatric disorders present",
        "Overweight",
        "Parkinson Disease",
        "Parkinson's",
        "Parkinson's disease",
        "Parkinsonism",
        "Parkinsonism clinical diagnosis",
        "Parkinsonism with asymmetrical onset",
        "Parkinson’s disease",
        "Pelvic Dystonia Present",
        "Pelvis tremor present",
        "Pelvis tremor type",
        "Polycythemia",
        "Polyps",
        "Postural Instability",
        "Primary Dystonia type",
        "Primary clinical diagnosis",
        "Progressive supranuclear palsy",
        "Prolonged febrile seizure",
        "Psychiatric disorder present",
        "Resting Tremor",
        "Rheumatoid arthritis",
        "Rigidity",
        "Schizophrenia",
        "Scoliosis",
        "Secondary neurological diagnoses",
        "Segmental Dystonia",
        "Seizure Disorders - Acute Symptomatic Seizures (Age of onset, if applicable)",
        "Seizure Disorders - Epilepsy: >2 unprovoked seizures (Age of onset, if applicable)",
        "Seizure Disorders - Febrile seizures (Age of onset, if applicable)",
        "Seizure Disorders - Other (Age of onset, if applicable) (Age of onset, if applicable)",
        "Seizure Disorders - Single Unprovoked seizure (Age at occurrence, if applicable)",
        "Seizure Disorders - Status Epilepticus (Age of onset, if applicable)",
        "Severe weight loss",
        "Shortness of breath",
        "Shoulder Dystonia Present",
        "Shoulder tremor present",
        "Shoulder tremor type",
        "Significant cognitive impairment or dementia",
        "Simple febrile seizure",
        "Status epilepticus",
        "Stroke",
        "Symptoms: Motor Tics",
        "Symptoms: Vocal/Verbal Tics",
        "Thyroid (hypothyroidism)",
        "Tongue Dystonia present",
        "Tongue tremor",
        "Tongue tremor type",
        "Tonic",
        "Tonic-Clonic",
        "Tourette Syndrome or tic disorders",
        "Tourettes",
        "Tourette’s Disorder",
        "Tremor",
        "Trunk Dystonia present",
        "Trunk tremor present",
        "Trunk tremor type",
        "Ulcers",
        "Unclassifiable seizure type",
        "Upper Arm Dystonia present",
        "Upper Arm tremor present",
        "Upper Arm tremor type",
        "Upper Face Dystonia present",
        "Upper Face tremor present",
        "Upper Face tremor type",
        "Upper Leg Dystonia present",
        "Upper Leg tremor present",
        "Upper Leg tremor type",
        "Urinary frequency",
        "Valvular disease",
        "Vasculitis",
        "Vertical gaze palsy",
        "Whipple's disease",
        "atrial fibrillation",
        "axial rigidity",
        "diabetes mellitus",
        "diagnosis",
        "dysautonomia",
        "gaze palsy",
        "hallucinations",
        "heart disease",
        "hypertension",
        "ischemic stroke subtype based on TOAST",
        "memory loss",
        "myocardial infarction",
        "stroke",
    }
)

# Per the mapping: "associatedPhenotype.name" / "associatedPhenotype" columns.
# Values that are simply category labels (e.g. Race=White) are emitted as
# {"name": <value>} DefinedTerm objects.
_DEFINED_TERM_PHENOTYPE_FIELDS = (
    "Race",
    "Ethnicity",
    "Hispanic or Latino/Not Hispanic or Latino",
)

_ANATOMY_FIELDS = (
    "Sample Source",
    "Source",
    "Biopsy Source",
    "Tissue Type",
)

_LOCATION_FIELDS = (
    "Country of Origin",
    "Country",
)


_KEYWORDS_FIELDS = ("Subcollection",)


_QUANTITY_RE = re.compile(r"^\s*(?P<value>[-+]?\d+(?:\.\d+)?)\s*(?P<unit>.+?)\s*$")


def _parse_quantity(raw):
    """Parse '10 µg' -> (10.0, 'µg'); '25 YR' -> (25.0, 'YR'); else (None, raw)."""
    if not raw:
        return None, None
    m = _QUANTITY_RE.match(str(raw))
    if not m:
        return None, str(raw).strip()
    try:
        value = float(m.group("value"))
    except (TypeError, ValueError):
        return None, str(raw).strip()
    unit = m.group("unit").strip() or None
    if unit and unit.casefold() in {"µg", "ug", "mcg"}:
        unit = "ug"
    return value, unit


_GEO_RE = re.compile(r"\b(GSM\d+|GSE\d+|GDS\d+|GPL\d+)\b", re.IGNORECASE)


def _extract_geo_id(raw):
    if not raw:
        return None
    m = _GEO_RE.search(str(raw))
    if m:
        return m.group(1).upper()
    return None


def _add_health_condition(output, name):
    name = _normalize(name)
    if not name:
        return
    insert_value(output, "healthCondition", {"name": name})


def _add_associated_genotype(output, value):
    value = _normalize(value)
    if not value:
        return
    insert_value(output, "associatedGenotype", value)


class USIDNETItemProcessorPipeline:
    def process_item(self, item, spider=None):
        if not item:
            return None

        catalog_id = _normalize(item.get("CatalogID"))
        if not catalog_id:
            return None

        url = f"https://www.coriell.org/0/Sections/Search/Sample_Detail.aspx?Ref={catalog_id}"

        output = {
            "@context": "http://schema.org/",
            "@type": "Sample",
            "_id": f"{catalog_id}",
            "identifier": catalog_id,
            "url": url,
            "distribution": [{"@type": "DataDownload", "contentUrl": url}],
            "includedInDataCatalog": {
                "@type": "DataCatalog",
                "name": "USIDNET",
                "url": "https://www.usidnet.org/",
                "versionDate": datetime.date.today().isoformat(),
                "archivedAt": url,
            },
            "additionalType": "BioSample",
            "creativeWorkStatus": "Available",
            "conditionsOfAccess": "Restricted",
            "sampleAvailability": True,
        }

        # ---- isAccessibleForFree (Coriell shows tiered pricing; any $0.00 tier => free)
        amounts = []
        for price in item.get("Prices") or []:
            m = re.search(r"[-+]?\d+(?:\.\d+)?", str(price).replace(",", ""))
            if m:
                amounts.append(float(m.group()))
        if amounts:
            output["isAccessibleForFree"] = any(a == 0 for a in amounts)

        # ---- sampleType (prefer the descriptive "Product", fall back to "ProductTypeID")
        if product := _normalize(item.get("Product")):
            insert_value(output, "sampleType", {"name": product})
        elif product_type_id := _normalize(item.get("ProductTypeID")):
            insert_value(output, "sampleType", {"name": product_type_id})

        # ---- sdPublisher
        if sdpublisher := _normalize(item.get("Repository")):
            insert_value(output, "sdPublisher", {"@type": "Organization", "name": sdpublisher})

        # ---- description (single string; concat Description + Remarks)
        for field in ("Description", "Remarks"):
            val = _normalize(item.get(field))
            if val:
                insert_value(output, "description", val, extend=True)

        # ---- species / infectiousAgent (mapping says both)
        species_name = _normalize(item.get("Species"))
        common_name = _normalize(item.get("Common Name"))
        if species_name:
            sp = {"name": species_name}
            if common_name:
                sp["commonName"] = common_name
            insert_value(output, "species", sp)
            ia = {"name": species_name}
            if common_name:
                ia["commonName"] = common_name
            insert_value(output, "infectiousAgent", ia)

        # ---- sex
        sex = _normalize(item.get("Sex"))
        if sex:
            insert_value(output, "sex", sex)

        # ---- keywords
        if keywords := _normalize(item.get("Subcollection")):
            insert_value(output, "keywords", keywords)

        # ---- anatomicalStructure
        for field in _ANATOMY_FIELDS:
            val = _normalize(item.get(field))
            if val and not _is_excluded(val):
                insert_value(output, "anatomicalStructure", {"name": val})

        # ---- developmentalStage
        # "Age at Sampling" takes precedence; otherwise "Age".
        # "Age (Yrs)" is a numeric column → QuantitativeValue with unitText=years.
        age_at_sampling = _normalize(item.get("Age at Sampling"))
        age = _normalize(item.get("Age"))
        age_yrs = _normalize(item.get("Age (Yrs)"))
        dev_stage: dict = {}
        if age_at_sampling and not _is_excluded(age_at_sampling):
            value, unit = _parse_quantity(age_at_sampling)
            if value is not None:
                dev_stage["value"] = value
                if unit:
                    dev_stage["unitText"] = unit
            dev_stage["name"] = age_at_sampling
        elif age and not _is_excluded(age):
            dev_stage["name"] = age
        if age_yrs and not _is_excluded(age_yrs):
            try:
                dev_stage["value"] = float(age_yrs)
                dev_stage.setdefault("unitText", "years")
                dev_stage.setdefault("name", f"{age_yrs} years")
            except (TypeError, ValueError):
                logger.warning(f"Could not parse age in years: {age_yrs}")
        if dev_stage:
            insert_value(output, "developmentalStage", dev_stage)

        # ---- associatedPhenotype (DefinedTerm-style fields)
        for field in _DEFINED_TERM_PHENOTYPE_FIELDS:
            val = _normalize(item.get(field))
            if val and not _is_excluded(val):
                insert_value(output, "associatedPhenotype", {"name": val})

        # ---- associatedPhenotype: Handedness ("Right" -> "Right Handedness")
        for field in ("Handedness", "Handedness:"):
            val = _normalize(item.get(field))
            if val and not _is_excluded(val):
                # Spider may concatenate multiple radio options into a single string
                # ("Right   Left  Ambidextrous"). Best-effort: pick the first token.
                first = val.split()[0]
                insert_value(output, "associatedPhenotype", {"name": f"{first} Handedness"})

        # ---- associatedPhenotype: Unilateral Babinski sign  ("{property}{value}")
        ubs = _normalize(item.get("Unilateral Babinski sign"))
        if ubs and not _is_excluded(ubs):
            insert_value(output, "associatedPhenotype", {"name": f"Unilateral Babinski sign {ubs}"})

        # ---- sampleQuantity ("Quantity": "10 µg")
        quantity = _normalize(item.get("Quantity"))
        if quantity and not _is_excluded(quantity):
            value, unit = _parse_quantity(quantity)
            sq: dict = {"name": quantity}
            if value is not None:
                sq["value"] = value
            if unit:
                sq["unitText"] = unit
            insert_value(output, "sampleQuantity", sq)

        # ---- locationOfOrigin (AdministrativeArea)
        for field in _LOCATION_FIELDS:
            val = _normalize(item.get(field))
            if val and not _is_excluded(val):
                insert_value(
                    output,
                    "locationOfOrigin",
                    {"@type": "AdministrativeArea", "name": val},
                )

        # ---- cellType (DefinedTerm) — accept "Cell Type" and "Cell Subtype"
        for field in ("Cell Type", "Cell Subtype"):
            val = _normalize(item.get(field))
            if val and not _is_excluded(val):
                insert_value(output, "cellType", {"name": val})

        # ---- sampleProcess
        ipsc = _normalize(item.get("Induced Pluripotent Stem Cell"))
        if ipsc and not _is_excluded(ipsc):
            insert_value(output, "sampleProcess", ipsc)

        # ---- sampleState (text)
        passage = _normalize(item.get("Passage Frozen"))
        if passage and not _is_excluded(passage):
            insert_value(output, "sampleState", f"{passage} passage")
        pdl = _normalize(item.get("PDL at Freeze"))
        if pdl and not _is_excluded(pdl):
            insert_value(output, "sampleState", f"{pdl} passage doubling level")

        # ---- associatedGenotype
        gene = _normalize(item.get("Gene"))
        chrom_loc = _normalize(item.get("Chromosomal Location"))
        mutations = _normalize(item.get("Mutations"))
        identified_mutation = _normalize(item.get("Identified Mutation"))
        if gene or chrom_loc or mutations:
            primary = " ".join(p for p in (gene, chrom_loc, mutations) if p).strip()
            if primary:
                _add_associated_genotype(output, primary)
        if (
            identified_mutation
            and identified_mutation != mutations
            and not _is_excluded(identified_mutation)
        ):
            alt = " ".join(p for p in (gene, chrom_loc, identified_mutation) if p).strip()
            if alt:
                _add_associated_genotype(output, alt)
        iscn = _normalize(item.get("ISCN"))
        if iscn and not _is_excluded(iscn):
            _add_associated_genotype(output, iscn)
        fragile_x = _normalize(item.get("Fragile X"))
        if fragile_x and not _is_excluded(fragile_x):
            _add_associated_genotype(output, f"Fragile X: {fragile_x}")

        # ---- isBasisFor.identifier (GEO accession)
        geo = _normalize(item.get("GEO"))
        geo_id = _extract_geo_id(geo)
        if geo_id:
            insert_value(output, "isBasisFor", {"identifier": geo_id})

        # ---- healthCondition (long-tail of disease columns from the mapping)
        for hc_field in _HEALTH_CONDITION_FIELDS:
            if hc_field not in item:
                continue
            raw = item.get(hc_field)
            if _is_excluded(raw):
                continue
            value = _normalize(raw)
            if hc_field in _HC_VALUE_IS_DISEASE:
                _add_health_condition(output, value)
            else:
                _add_health_condition(output, hc_field)

        return output
