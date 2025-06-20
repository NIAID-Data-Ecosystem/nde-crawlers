import json
import os

import text2term


def read_ndjson(file_path):
    docs = []
    with open(file_path, "r") as file:
        for line in file:
            docs.append(json.loads(line.strip()))
    return docs


def add_topic_category(docs, source_name):
    """
    Adds 'topicCategory' to documents and appends new topics if they are not duplicates
    (based on the 'identifier' property). Loads topic categories from a specified JSON file
    in the /data/nde-hub/topic_categories/ directory.
    """
    if isinstance(docs, str):
        file_path = os.path.join(docs, "data.ndjson")
        docs = read_ndjson(file_path)

    # Cache the EDAM ontology if not already cached
    if not os.path.exists("/data/nde-hub/topic_categories/cache/edam"):
        text2term.cache_ontology("https://edamontology.org/EDAM_unstable.owl", "edam")

    file_path = f"/data/nde-hub/topic_categories/{source_name}.json"
    with open(file_path, "r") as file:
        topic_categories = json.load(file)

    topic_dict = {list(item.keys())[0]: list(item.values())[0] for item in topic_categories}

    all_topics = set()
    for topics in topic_dict.values():
        all_topics.update([topic.strip('"') for topic in topics])
    all_topics = list(all_topics)

    # Batch query the topics
    t2t_result = text2term.map_terms(all_topics, "edam", use_cache=True, base_iris=["http://edamontology.org/topic_"])
    t2t_result.sort_values(["Source Term", "Mapping Score"], ascending=[True, False], inplace=True)

    # Create a mapping for exact matches only
    topic_mapping = {}
    for _, row in t2t_result.iterrows():
        source_term = row["Source Term"]
        mapped_term_label = row["Mapped Term Label"]
        if source_term.lower() == mapped_term_label.lower():
            topic_mapping[source_term] = row

    # Define the exception topics
    exception_topics = {"Anatomy", "Transcriptomics", "Developmental biology", "Oncology", "Physiology"}

    updated_docs = []
    for doc in docs:
        doc_id = doc["_id"].lower()
        topics = topic_dict.get(doc_id, [])
        if topics:
            # Ensure topicCategory exists; if not, initialize it.
            if "topicCategory" not in doc:
                doc["topicCategory"] = []
            # Create a set of existing identifiers for de-duplication.
            existing_ids = {topic_item.get("identifier") for topic_item in doc["topicCategory"]}

            contains_exception_topic = any(topic.strip('"') in exception_topics for topic in topics)
            for topic in topics:
                topic_cleaned = topic.strip('"')
                if topic_cleaned in topic_mapping:
                    if topic_cleaned == "Human biology" and not contains_exception_topic:
                        continue  # Skip 'Human biology' if no exception topic is present
                    row = topic_mapping[topic_cleaned]
                    identifier = "topic_" + row["Mapped Term CURIE"].split(":")[1]
                    # Only add if this identifier is not already present
                    if identifier not in existing_ids:
                        new_topic = {
                            "@type": "DefinedTerm",
                            "name": topic_cleaned,
                            "curatedBy": {"name": "GPT-4o-mini", "url": "https://openai.com/index/chatgpt"},
                            "identifier": identifier,
                            "url": row["Mapped Term IRI"],
                            "inDefinedTermSet": "EDAM",
                            "fromGPT": True,
                        }
                        doc["topicCategory"].append(new_topic)
                        existing_ids.add(identifier)
        updated_docs.append(doc)
    return updated_docs
