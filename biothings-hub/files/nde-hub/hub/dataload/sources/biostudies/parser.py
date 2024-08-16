import json
import os

import dateutil


def parse_file(doc):
    accno = doc.pop("accno")

    output = {"_id": f"biostudies_{accno}", "url": f"https://www.ebi.ac.uk/studies/{accno}"}

    # handle attributes
    if attributes := doc.pop("attributes", None):
        for attribute in attributes:
            key = attribute.pop("name").casefold()
            value = attribute.pop("value")

            known_attributes = ["template", "rootpath"]

            if key == "attachto":
                output["url"] = f"https://www.ebi.ac.uk/{value}/studies/{accno}"
            elif key == "releasedate":
                output["releaseDate"] = dateutil.parser.parse(value, ignoretz=True).date().isoformat()
            elif key == "doi":
                output["doi"] = value
            elif key == "title":
                output["name"] = value
            elif key in known_attributes:
                pass
            else:
                print(f"Accno: {accno}, Attribute name: {key}")

    return output


def parse_file_dir(input_dir):
    for filename in os.listdir(input_dir):
        if filename.endswith(".json"):
            with open(os.path.join(input_dir, filename), "r") as infile:
                doc = json.load(infile)
                yield parse_file(doc)
