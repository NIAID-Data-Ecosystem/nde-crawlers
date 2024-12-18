import csv
import requests
import os
from pymongo import MongoClient
from datetime import datetime

# Configuration
MONGO_URI = os.environ.get("MONGO_URI")
DATABASE_NAME = "nde_hub"
OUTPUT_CSV = "ncbi_linkout.csv"

PROVIDER_ID = "4060"
DB_NAME = "PubMed"
BASE_URL = "https://data.niaid.nih.gov/resources?id="
ICON_URL = ""
URL_NAME = "NIAID Data Ecosystem"
SUBJECT_TYPE = ""
ATTRIBUTE = ""

def get_latest_metadata():
    r = requests.get("https://api.data.niaid.nih.gov/v1/metadata")
    r.raise_for_status()
    return r.json()

def get_pmid_documents(collection, size):
    return collection.find({
        "citation": {
            "$elemMatch": {
                "pmid": {"$exists": True, "$ne": None}
            }
        }
    }).limit(size)


def main():
    metadata = get_latest_metadata()
    build_date = metadata.get("build_date")
    build_version = metadata.get("build_version")

    # Check last build date if needed
    if os.path.exists("last_build_date.txt"):
        with open("last_build_date.txt", "r") as f:
            last_build_date = f.read().strip()
    else:
        last_build_date = None

    if last_build_date == build_date:
        print("No new build. Exiting.")
        return

    client = MongoClient(MONGO_URI)
    db = client[DATABASE_NAME]

    # Find any collection that starts with nde_all_prod_<build_version>
    pattern = f"nde_all_prod_{build_version}"
    all_collections = db.list_collection_names()
    target_collection = None
    for coll in all_collections:
        if coll.startswith(pattern):
            target_collection = coll
            break

    if not target_collection:
        print("Could not find the target collection.")
        return

    # Generate the CSV
    size = 50
    with open(OUTPUT_CSV, "w", newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["PrId", "DB", "UID", "URL", "IconUrl", "UrlName", "SubjectType", "Attribute"])

        for doc in get_pmid_documents(db[target_collection], size):
            citations = doc.get("citation", [])
            dataset_id = str(doc.get("_id", ""))
            for c in citations:
                pmid = c.get("pmid")
                if pmid:
                    url = f"{BASE_URL}{dataset_id}"
                    writer.writerow([PROVIDER_ID, DB_NAME, pmid, url, ICON_URL, URL_NAME, SUBJECT_TYPE, ATTRIBUTE])

    print(f"CSV file generated: {OUTPUT_CSV}")

    # Update last_build_date.txt
    with open("last_build_date.txt", "w") as f:
        f.write(build_date)

if __name__ == "__main__":
    main()
