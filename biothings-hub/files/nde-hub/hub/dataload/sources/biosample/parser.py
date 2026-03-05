import json
import os
from pathlib import Path

import orjson
from utils.sex import add_sex


def parse_sex(data_folder):
    """
    Parse a GEO SOFT platform file into a dictionary.
    Each key is the SOFT field (e.g., '!Platform_title'), value is a list if repeated, or a string.
    """
    with open(Path(__file__).resolve().parent / "sex_mappings.json", "r") as f:
        sex_mapping = json.load(f)

    with open(os.path.join(data_folder, "data.ndjson"), "rb") as f:
        for line in f:
            doc = orjson.loads(line)
            sex = doc.pop("sex", None)
            if sex:
                if not isinstance(sex, list):
                    sex = [sex]
                for s in sex:
                    add_sex(s, doc, sex_mapping)
            yield doc

