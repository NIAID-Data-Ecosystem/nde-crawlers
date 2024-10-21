import os
import sqlite3

import orjson
import text2term
from config import logger

DB_PATH = "/data/nde-hub/standardizers/measurement_technique_lookup/measurement_technique_lookup.db"
CACHE = {}  # In-memory cache for repeated lookups


def create_sqlite_db(conn):
    """
    Create the SQLite database tables if they don't exist.
    """
    with conn:
        conn.execute(
            """CREATE TABLE IF NOT EXISTS measurement_technique_lookup (
                 technique_name TEXT PRIMARY KEY,
                 standardized_data TEXT
               )"""
        )


def update_sqlite_db_bulk(conn, batch):
    """
    Perform a bulk update on the SQLite database.
    """
    with conn:
        conn.executemany(
            "INSERT OR REPLACE INTO measurement_technique_lookup (technique_name, standardized_data) VALUES (?, ?)",
            batch,
        )


def sqlite_lookup(conn, technique_name):
    """
    Look up the standardized data for a given measurement technique.
    Uses in-memory cache for faster repeated lookups.
    """
    if technique_name in CACHE:
        return CACHE[technique_name]

    cur = conn.execute(
        "SELECT standardized_data FROM measurement_technique_lookup WHERE technique_name=?", (technique_name.lower(),)
    )
    result = cur.fetchone()
    if result:
        data = orjson.loads(result[0])
        CACHE[technique_name] = data
        return data
    return None


def find_best_ontology_match(technique_name):
    """
    Use text2term to find the best ontology match for the measurement technique.
    Caches results to avoid redundant lookups.
    """
    best_match = None
    highest_score = 0

    ontologies = ["edam", "mmo", "ncit", "efo", "bao", "chmo", "obi"]

    for ontology in ontologies:
        try:
            matches = text2term.map_terms([technique_name], ontology, use_cache=True)
        except Exception as e:
            logger.error(f"Error mapping technique to ontology '{ontology}': {e}")
            continue

        if not matches.empty:
            top_match = matches.sort_values(by="Mapping Score", ascending=False).iloc[0]
            score = top_match["Mapping Score"]
            if score > highest_score:
                highest_score = score
                iri = top_match["Mapped Term IRI"]
                curie = top_match["Mapped Term CURIE"]
                best_match = {
                    "name": top_match["Mapped Term Label"],
                    "originalName": technique_name,
                    "identifier": curie.split(":")[1],
                    "inDefinedTermSet": curie.split(":")[0],
                    "isCurated": True,
                    "url": iri,
                    "curie": curie,
                    "score": float(score),
                }

    if best_match:
        return best_match
    else:
        return {
            "name": technique_name,
        }


def process_measurement_technique(data):
    """
    Process and standardize measurement techniques in the document list.
    Optimized to use in-memory caching and batch database updates.
    """
    conn = sqlite3.connect(DB_PATH)
    create_sqlite_db(conn)

    if isinstance(data, str):
        with open(os.path.join(data, "data.ndjson"), "rb") as f:
            doc_list = [orjson.loads(line) for line in f]
    else:
        doc_list = list(data)

    count = 0
    update_batch = []  # Batch for bulk inserts
    for doc in doc_list:
        count += 1
        if count % 1000 == 0:
            logger.info(f"Processing document {count}...")

        if "measurementTechnique" in doc:
            measurement_data = doc["measurementTechnique"]
            techniques = [measurement_data] if isinstance(measurement_data, dict) else measurement_data

            standardized_techniques = []
            for technique_obj in techniques:
                technique_name = technique_obj.get("name", "")
                if not technique_name:
                    logger.warning(f"No name found in measurementTechnique object: {technique_obj}")
                    continue

                standardized_data = sqlite_lookup(conn, technique_name)
                if not standardized_data:
                    standardized_data = find_best_ontology_match(technique_name)
                    update_batch.append((technique_name.lower(), orjson.dumps(standardized_data).decode("utf-8")))
                    CACHE[technique_name] = standardized_data

                standardized_techniques.append(standardized_data)

            doc["measurementTechnique"] = (
                standardized_techniques[0] if len(standardized_techniques) == 1 else standardized_techniques
            )

        if len(update_batch) >= 500:
            update_sqlite_db_bulk(conn, update_batch)
            update_batch = []

    if update_batch:
        update_sqlite_db_bulk(conn, update_batch)

    conn.close()
    return doc_list
