from hub.dataload.nde import NDESourceUploader

COVID_HEALTH_CONDITION = {
    "@type": "DefinedTerm",
    "alternateName": [
        "2019 novel coronavirus infection",
        "2019-nCoV infection",
        "coronavirus disease 2019",
        "severe acute respiratory syndrome coronavirus 2 infectious disease",
    ],
    "curatedBy": {
        "@type": "SoftwareApplication",
        "dateModified": "2025-01-15",
        "name": "Biothings API",
        "url": "https://biothings.io/",
    },
    "identifier": "0100096",
    "inDefinedTermSet": "MONDO",
    "isCurated": True,
    "name": "COVID-19",
    "originalName": "COVID-19",
    "url": "http://purl.obolibrary.org/obo/MONDO_0100096",
}

HOMO_SAPIENS = {
    "@type": "DefinedTerm",
    "identifier": "9606",
    "inDefinedTermSet": "UniProt",
    "url": "https://www.uniprot.org/taxonomy/9606",
    "originalName": "homo sapiens",
    "isCurated": True,
    "curatedBy": {
        "@type": "SoftwareApplication",
        "name": "PubTator",
        "url": "https://www.ncbi.nlm.nih.gov/research/pubtator/api.html",
        "dateModified": "2023-10-05",
    },
    "name": "Homo sapiens",
    "commonName": "Human",
    "displayName": "Human | Homo sapiens",
    "alternateName": [
        "Human",
        "Homo sapiens Linnaeus, 1758",
        "human",
        "Home sapiens",
        "Homo sampiens",
        "Homo sapeins",
        "Homo sapian",
        "Homo sapians",
        "Homo sapien",
        "Homo sapience",
        "Homo sapiense",
        "Homo sapients",
        "Homo sapines",
        "Homo spaiens",
        "Homo spiens",
        "Humo sapiens",
    ],
    "classification": "host",
}

COVID_INFECTIOUS_AGENT = {
    "@type": "DefinedTerm",
    "alternateName": [
        "2019-nCoV",
        "Wuhan coronavirus",
        "SARS-2",
        "SARS-CoV2",
        "Wuhan seafood market pneumonia virus",
        "HCoV-19",
        "COVID19",
        "COVID-19 virus",
        "Human coronavirus 2019",
        "COVID-19",
    ],
    "classification": "infectiousAgent",
    "commonName": "2019-nCoV",
    "curatedBy": {
        "@type": "SoftwareApplication",
        "dateModified": "2025-02-09",
        "name": "Data Discovery Engine",
        "url": "https://discovery.biothings.io/",
    },
    "displayName": "2019-nCoV | Severe acute respiratory syndrome coronavirus 2",
    "identifier": "2697049",
    "inDefinedTermSet": "UniProt",
    "isCurated": True,
    "name": "Severe acute respiratory syndrome coronavirus 2",
    "originalName": "Severe acute respiratory syndrome coronavirus 2",
    "url": "https://www.uniprot.org/taxonomy/2697049",
}


class Covid_Radx_Uploader(NDESourceUploader):
    name = "covid_radx"

    def post_process(self, doc):
        """Every RADx record is COVID-19 in humans, whatever the source metadata says."""
        # Make sure healthCondition and infectiousAgent are lists.
        for field in ["healthCondition", "infectiousAgent", "species"]:
            if field not in doc:
                doc[field] = []
            elif not isinstance(doc[field], list):
                doc[field] = [doc[field]]

        # Check for existing COVID-19 health condition (by identifier)
        if not any(isinstance(item, dict) and item.get("identifier") == "0100096" for item in doc["healthCondition"]):
            doc["healthCondition"].append(dict(COVID_HEALTH_CONDITION))

        # Check for existing SARS-CoV-2 infectious agent (by identifier)
        if not any(isinstance(item, dict) and item.get("identifier") == "2697049" for item in doc["infectiousAgent"]):
            doc["infectiousAgent"].append(dict(COVID_INFECTIOUS_AGENT))

        # Overwrite the species field to homo sapiens
        doc["species"] = [dict(HOMO_SAPIENS)]
        return doc
