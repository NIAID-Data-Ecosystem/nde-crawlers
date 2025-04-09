from hub.dataload.nde import NDESourceUploader
from utils.corrections import corrections
from utils.extract import process_descriptions
from utils.funding_helper import standardize_funding
from utils.measurement_technique_helper import process_measurement_technique
from utils.pmid_helper import load_pmid_ctfd
from utils.pubtator import standardize_data
from utils.topic_category_helper import add_topic_category
from utils.utils import nde_upload_wrapper


class Covid_Radx_Uploader(NDESourceUploader):
    name = "covid_radx"

    @nde_upload_wrapper
    def load_data(self, data_folder):
        docs = load_pmid_ctfd(data_folder)
        docs = standardize_funding(docs)
        docs = standardize_data(docs)
        docs = process_descriptions(docs)
        docs = corrections(docs)
        docs = process_measurement_technique(docs, self.name)
        docs = add_topic_category(docs, self.name)

        covid_health_condition = {
            "alternateName": [
                "2019 novel coronavirus infection",
                "2019-nCoV infection",
                "coronavirus disease 2019",
                "severe acute respiratory syndrome coronavirus 2 infectious disease",
            ],
            "curatedBy": {"dateModified": "2025-01-15", "name": "Biothings API", "url": "https://biothings.io/"},
            "identifier": "0100096",
            "inDefinedTermSet": "MONDO",
            "isCurated": True,
            "name": "COVID-19",
            "originalName": "COVID-19",
            "url": "http://purl.obolibrary.org/obo/MONDO_0100096",
        }

        homo_sapiens = {
            "identifier": "9606",
            "inDefinedTermSet": "UniProt",
            "url": "https://www.uniprot.org/taxonomy/9606",
            "originalName": "homo sapiens",
            "isCurated": True,
            "curatedBy": {
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

        covid_infectious_agent = {
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

        for doc in docs:
            # Make sure healthCondition and infectiousAgent are lists.
            for field in ["healthCondition", "infectiousAgent", "species"]:
                if field not in doc:
                    doc[field] = []
                elif not isinstance(doc[field], list):
                    doc[field] = [doc[field]]

            # Check for existing COVID-19 health condition (by identifier)
            if not any(
                isinstance(item, dict) and item.get("identifier") == "0100096" for item in doc["healthCondition"]
            ):
                doc["healthCondition"].append(covid_health_condition)

            # Check for existing SARS-CoV-2 infectious agent (by identifier)
            if not any(
                isinstance(item, dict) and item.get("identifier") == "2697049" for item in doc["infectiousAgent"]
            ):
                doc["infectiousAgent"].append(covid_infectious_agent)

            # Overwrite the species field to homo sapiens
            doc["species"] = [homo_sapiens]

            yield doc
