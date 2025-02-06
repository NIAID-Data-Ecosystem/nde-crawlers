COMPUTATIONAL_TOOL_REQUIRED = [
    "date",
    "includedInDataCatalog",
    "funding",
    "author",
    "description",
    "name"
]
COMPUTATIONAL_TOOL_RECOMMENDED = [
    "citedBy",
    "doi",
    "topicCategory",
    "codeRepository",
    "programmingLanguage",
    "applicationCategory",
    "applicationSubCategory",
    "input",
    "output",
    "featureList",
    "operatingSystem",
    "softwareRequirements",
    "softwareVersion",
    "citation",
    "conditionsOfAccess",
    "dateModified",
    "interactionStatistic",
    "license",
    "identifier",
    "url"
]
COMPUTATIONAL_TOOL_REQUIRED_AUGMENTED = ['funding']
COMPUTATIONAL_TOOL_RECOMMENDED_AUGMENTED = ["citation", "topicCategory"]
RESOURCE_CATALOG_REQUIRED = [
    "date",
    "funding",
    "includedInDataCatalog",
    "measurementTechnique",
    "description",
    "name",
    "url",
    "about",
    "genre",
    "collectionType"
]
RESOURCE_CATALOG_RECOMMENDED = [
    "author",
    "citedBy",
    "doi",
    "infectiousAgent",
    "healthCondition",
    "species",
    "variableMeasured",
    "citation",
    "conditionsOfAccess",
    "dateCreated",
    "dateModified",
    "datePublished",
    "interactionStatistic",
    "isBasedOn",
    "keywords",
    "license",
    "sdPublisher",
    "spatialCoverage",
    "temporalCoverage",
    "usageInfo",
    "identifier",
    "topicCategory",
    "collectionSize",
    "hasAPI",
    "hasDownload",
]
RESOURCE_CATALOG_REQUIRED_AUGMENTED = ["funding", "measurementTechnique"]
RESOURCE_CATALOG_RECOMMENDED_AUGMENTED = [
    "species",
    "infectiousAgent",
    "healthCondition",
    "citation",
    "topicCategory"
]
DATASET_REQUIRED_FIELDS = [
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
DATASET_RECOMMENDED_FIELDS = [
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
DATASET_REQUIRED_AUGMENTED_FIELDS = ["funding", "measurementTechnique"]
DATASET_RECOMMENDED_AUGMENTED_FIELDS = [
    "species", "infectiousAgent", "healthCondition", "citation", "topicCategory"]
