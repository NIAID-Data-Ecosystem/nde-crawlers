MAPPING_SCORES = {
    "abstract": 0.4,
    "alternateName": 0.4,
    "author": {
        "familyName": 0.2,
        "givenName": 0.2,
        "name": 0.4,
    },
    "citation": {
        "author": {
            "familyName": 0.1,
            "givenName": 0.1,
            "name": 0.2,
        },
        "doi": 0.2,
        "name": 0.1,
        "pmid": 0.1,
        "url": 0.1,
    },
    "citedBy": {
        "doi": 0.1,
        "name": 0.1,
        "pmid": 0.1,
        "url": 0.1,
    },
    "contentUrl": 0.2,
    "dateCreated": 0.2,
    "dateModified": 0.2,
    "datePublished": 0.2,
    "distribution": {
        "contentUrl": 0.1,
    },
    "doi": 0.3,
    "funding": {
        "url": 0.3,
        "name": 0.1,
        "funder": {
            "name": 0.1,
            "url": 0.3,
        },
    },
    "healthCondition": {
        "name": 0.3,
        "isCurated": 0.3,
    },
    "infectiousAgent": {
        "name": 0.3,
        "isCurated": 0.3,
    },
    "isBasedOn": {
        "doi": 0.3,
        "name": 0.1,
        "pmid": 0.3,
        "url": 0.1,
    },
    "keywords": 0.1,
    "measurementTechnique": {
        "name": 0.1,
    },
    "sdPublisher": {
        "name": 0.3,
        "url": 0.3,
    },
    "species": {
        "name": 0.3,
        "isCurated": 0.3,
    },
}
REQUIRED_FIELDS = [
    "name",
    "description",
    "author",
    "url",
    "measurementTechnique",
    "includedInDataCatalog",
    "distribution",
    "funding",
    "date"
]
RECOMMENDED_FIELDS = [
    "dateCreated",
    "dateModified",
    "datePublished",
    "citedBy",
    "doi",
    "infectiousAgent",
    "healthCondition",
    "species",
    "variableMeasured",
    "citation",
    "conditionsOfAccess",
    "isBasedOn",
    "keywords",
    "license",
    "sdPublisher",
    "spatialCoverage",
    "temporalCoverage",
    "topicCategory",
    "identifier",
    "usageInfo",
    "interactionStatistic"
]
REQUIRED_AUGMENTED_FIELDS = [
    "funding",
    "measurementTechnique"
]
RECOMMENDED_AUGMENTED_FIELDS = [
    "species",
    "infectiousAgent",
    "healthCondition",
    "citation"
]
