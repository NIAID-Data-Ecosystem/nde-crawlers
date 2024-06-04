import datetime
import ftplib
import json
import urllib.parse

import requests


def fetch_dataset_ids():
    base_url = "https://massive.ucsd.edu/ProteoSAFe/QueryDatasets"
    query_params = {
        "pageSize": 1000,
        "offset": 0,
        "query": urllib.parse.quote_plus('#{"query":{},"table_sort_history":"createdMillis_dsc"}'),
    }

    dataset_ids = []
    page = 1

    while True:
        response = requests.get(base_url, params=query_params)
        if response.status_code == 200:
            print(f"Processed page {page}")
            data = response.json()
            row_data = data.get("row_data", [])
            if not row_data:
                break  # No more data to process
            for item in row_data:
                dataset_ids.append(item["dataset"])
            query_params["offset"] += query_params["pageSize"]  # Increment offset for next page
            page += 1
            break
        else:
            raise Exception(f"Failed to fetch dataset IDs, status code: {response.status_code}")

    return dataset_ids

def fetch_dataset(dataset_id):
    url = f"https://massive.ucsd.edu/ProteoSAFe/QueryDatasets?pageSize=30&offset=0&query=%7B%22title_input%22:%22{dataset_id}%22%7D"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        raise Exception(f"Failed to fetch dataset {dataset_id}, status code: {response.status_code}")

def parse_dataset(json_data):
    output = {
        "includedInDataCatalog": {
            "name": "MassIVE",
            "url": "https://massive.ucsd.edu/ProteoSAFe/static/massive.jsp",
            "@type": "Dataset",
            "versionDate": datetime.datetime.today().strftime("%Y-%m-%d"),
        }
    }

    for item in json_data.get("row_data", []):
        if identifier := item.get("dataset"):
            output["identifier"] = identifier
            output['_id'] = "massive_" + identifier.lower()

        if task := item.get("task"):
            output["url"] = "https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=" + task

        if repo_path := item.get("repo_path"):
            path = repo_path.split("-")[-1]
            if path != "/data/massive":
                output["distribution"] = {"contentUrl": f"ftp://massive.ucsd.edu/{path}/{identifier}"}
            else:
                output["distribution"] = {"contentUrl": f"ftp://massive.ucsd.edu/v08/{identifier}"}

        if title := item.get("title"):
            output["name"] = title

        if sdPublisher := item.get("site"):
            output["sdPublisher"] = {"name": sdPublisher}

        if description := item.get("description"):
            output["description"] = description

        if keywords := item.get("keywords"):
            if "###" in keywords:
                output["keywords"] = keywords.split('###')
            else:
                output["keywords"] = keywords.split(', ')

        if date_created := item.get("create_time"):
            output["dateCreated"] = date_created

        if measurement_techniques := item.get("instrument_resolved"):
            mt_list = measurement_techniques.split("###")
            mt_names = [{"name": sp for sp in mt_list}]
            output["measurementTechnique"] = mt_names

        if species := item.get("species_resolved"):
            species_list = species.split("###")
            species_names = [{"name": sp.split(" (NCBITaxon:")[0]} for sp in species_list]
            output["species"] = species_names

        if author_name := item.get("pis"):
            pis_list = []
            for pi in author_name:
                pi_data = {}
                if name := pi.get("name"):
                    pi_data["name"] = name
                if email := pi.get("email"):
                    pi_data["email"] = email
                if institution := pi.get("institution"):
                    pi_data["affiliation"] = institution
                pis_list.append(pi_data)
            output["author"] = pis_list

        if publications := item.get("publications"):
            pub_list = []
            for pub in publications:
                pub_data = {}
                if citation_id := pub.get("id"):
                    pub_data["identifier"] = citation_id
                if authors := pub.get("authors"):
                    pub_data["author"] = authors
                if title := pub.get("title"):
                    pub_data["name"] = title
                if doi := pub.get("citation"):
                    pub_data["doi"] = doi
                if description := pub.get("abstract"):
                    pub_data["description"] = description
                pub_list.append(pub_data)
            output["citation"] = pub_list

        if privacy := item.get("privacy"):
            if privacy == "Public":
                output["conditionsOfAccess"] = "Open"
            else:
                output["conditionsOfAccess"] = "Closed"

        if 'measurementTechnique' not in output:
            output['measurementTechnique'] = [{"name": "Mass Spectrometry", "url": "https://ontobee.org/ontology/MMO?iri=http://purl.obolibrary.org/obo/MMO_0000534", "identifier": "MMO_0000534"}]
        else:
            output['measurementTechnique'].append({"name": "Mass Spectrometry", "url": "https://ontobee.org/ontology/MMO?iri=http://purl.obolibrary.org/obo/MMO_0000534", "identifier": "MMO_0000534"})


    return output


def parse():
    dataset_ids = fetch_dataset_ids()
    print(f"Found {len(dataset_ids)} datasets")
    count = 0
    for dataset_id in dataset_ids:
        count += 1
        if count % 100 == 0:
            print(f"Processed {count} datasets")
            break
        try:
            dataset_json = fetch_dataset(dataset_id)
            parsed_dataset = parse_dataset(dataset_json)
            yield parsed_dataset
        except Exception as e:
            print(f"Error processing dataset {dataset_id}: {e}")

records = []
for x in parse():
    records.append(x)
with open('massive.json', 'w') as f:
    json.dump(records, f)


# data = requests.get("https://massive.ucsd.edu/ProteoSAFe/QueryDatasets?pageSize=30&offset=0&query={%22title_input%22:%22MSV000094860%22}").json()['row_data'][0]

# # Prepare the CSV rows
# rows = []

# def extract_values(obj, prefix=''):
#     if isinstance(obj, dict):
#         for key, value in obj.items():
#             if isinstance(value, list):
#                 if len(value) > 0:
#                     extract_values(value[0], f"{prefix}{key}.")
#             else:
#                 rows.append([f"{prefix}{key}", value])
#     elif isinstance(obj, list):
#         for item in obj:
#             extract_values(item, prefix)

# # Extract values from the JSON
# extract_values(data)

# # Write to CSV
# with open('output.csv', 'w', newline='') as csvfile:
#     writer = csv.writer(csvfile)
#     writer.writerow(['property', 'value'])
#     writer.writerows(rows)

# print("CSV file has been created.")
