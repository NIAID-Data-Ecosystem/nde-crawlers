import csv
import os

import requests
from pymongo import MongoClient

# Configuration
# MONGO_URI = os.environ.get("MONGO_URI")
MONGO_URI = "mongodb://su11:27017"
DATABASE_NAME = "nde_hub"

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

def get_pmid_documents(collection):
    return collection.find({
        "citation": {
            "$elemMatch": {
                "pmid": {"$exists": True, "$ne": None}
            }
        },
        "includedInDataCatalog.name": {
            "$ne": "Omics Discovery Index (OmicsDI)"
        }
    })


def main():
    metadata = get_latest_metadata()
    build_date = metadata.get("build_date")
    build_date_str = build_date.split("T")[0].replace("-", "")  # "20241224"
    OUTPUT_CSV = f"ncbi_linkout_{build_date_str}.csv"
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
    with open(OUTPUT_CSV, "w", newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["PrId", "DB", "UID", "URL", "IconUrl", "UrlName", "SubjectType", "Attribute"])

        for doc in get_pmid_documents(db[target_collection]):
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
