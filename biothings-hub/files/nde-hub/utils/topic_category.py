"""Curated topicCategory terms.

Runs when the source has a topic file in `/data/nde-hub/topic_categories/`,
which lists the GPT-assigned topics per record id. Topics are matched to EDAM
with text2term (exact label matches only) and added as DefinedTerms, skipping
identifiers the record already carries.
"""

import json
import os
from functools import cache

import text2term
from config import logger

LOOKUP_DIR = "/data/nde-hub/topic_categories"
CACHE_DIR = f"{LOOKUP_DIR}/cache"
EDAM_URL = "https://edamontology.org/EDAM_unstable.owl"
EDAM_TOPIC_IRI = "http://edamontology.org/topic_"

CURATED_BY = {"name": "GPT-4o-mini", "url": "https://openai.com/index/chatgpt"}

# "Human biology" is only meaningful alongside one of these.
EXCEPTION_TOPICS = frozenset({"Anatomy", "Transcriptomics", "Developmental biology", "Oncology", "Physiology"})


def lookup_file(source):
    return f"{LOOKUP_DIR}/{source}.json"


@cache
def load_topics(source):
    """Load a source's topics: ({record id: [topic]}, {topic: EDAM term}).

    Cached, because mapping the source's whole topic vocabulary to EDAM is a
    single text2term run we only want to pay for once per upload.
    """
    with open(lookup_file(source), "r") as file:
        topic_categories = json.load(file)

    topic_dict = {list(item.keys())[0]: list(item.values())[0] for item in topic_categories}
    all_topics = sorted({topic.strip('"') for topics in topic_dict.values() for topic in topics})

    if not os.path.exists(os.path.join(CACHE_DIR, "edam")):
        logger.info("Caching the EDAM ontology in %s", CACHE_DIR)
        text2term.cache_ontology(EDAM_URL, "edam", cache_folder=CACHE_DIR)

    result = text2term.map_terms(
        all_topics,
        "edam",
        use_cache=True,
        base_iris=[EDAM_TOPIC_IRI],
        cache_folder=CACHE_DIR,
    )
    result.sort_values(["Source Term", "Mapping Score"], ascending=[True, False], inplace=True)

    # Exact label matches only.
    topic_mapping = {}
    for _, row in result.iterrows():
        source_term = row["Source Term"]
        if source_term.lower() == row["Mapped Term Label"].lower() and source_term not in topic_mapping:
            topic_mapping[source_term] = {
                "identifier": "topic_" + row["Mapped Term CURIE"].split(":")[1],
                "url": row["Mapped Term IRI"],
            }

    logger.info("Loaded topics for %s: %s records, %s EDAM terms", source, len(topic_dict), len(topic_mapping))
    return topic_dict, topic_mapping


def add_topic_category(docs, source):
    """Add curated topicCategory terms to every record in one batch."""
    topic_dict, topic_mapping = load_topics(source)

    for doc in docs:
        topics = topic_dict.get(doc["_id"].lower())
        if not topics:
            yield doc
            continue

        topic_categories = doc.setdefault("topicCategory", [])
        existing_ids = {topic_item.get("identifier") for topic_item in topic_categories}
        has_exception_topic = any(topic.strip('"') in EXCEPTION_TOPICS for topic in topics)

        for topic in topics:
            topic = topic.strip('"')
            term = topic_mapping.get(topic)
            if not term:
                continue
            if topic == "Human biology" and not has_exception_topic:
                continue
            if term["identifier"] in existing_ids:
                continue
            topic_categories.append(
                {
                    "@type": "DefinedTerm",
                    "name": topic,
                    "curatedBy": CURATED_BY,
                    "identifier": term["identifier"],
                    "url": term["url"],
                    "inDefinedTermSet": "EDAM",
                    "fromGPT": True,
                }
            )
            existing_ids.add(term["identifier"])
        yield doc
