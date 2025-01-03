REQUIRED_FIELDS = [
    "name",
    "description",
    "author",
    "url",
    "measurementTechnique",
    "includedInDataCatalog",
    "distribution",
    "funding",
    "date",
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
    "interactionStatistic",
]
REQUIRED_AUGMENTED_FIELDS = ["funding", "measurementTechnique"]
RECOMMENDED_AUGMENTED_FIELDS = ["species", "infectiousAgent", "healthCondition", "citation", "topicCategory"]
