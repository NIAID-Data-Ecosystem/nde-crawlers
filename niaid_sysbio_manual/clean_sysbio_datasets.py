import pandas as pd
from datetime import date
import requests

# Import the data from the template file
sheet_id = "12fm2Qkh39Ben7QS5QtIElVlaYwYY8RUExaFX9EX6iAg"
sheet_name = "2022-05-05_sysbio_datasets"
url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"

df_raw = pd.read_csv(url)
df = df_raw.loc[df_raw["Checked by RS?"] == "yes", :]



# nest related properties into a single object, from the flat columns
def nestAuthor(row):
    people = str2list(row["author.name"])
    affiliations = str2list(row["author.affiliation"])
    # interleave the two values together
    if(people is None):
        raise Exception(f"Missing author name at row {row}")
    elif(affiliations is None):
        # No affiliations to report
        return([{"name": person} for person in people])
    elif(len(people) == len(affiliations)):
        # paired name + affiliations
        return([{"name": pair[0], "dateModified": pair[1]} for pair in zip(people, affiliations)])
    elif(len(affiliations) == 1):
        # Assume affiliation applies to all names
        return([{"name": person, "affiliation": affiliations[0]} for person in people])
    else:
        # Attach array of affiliations to the person
        return([{"name": people[0], "affiliation": affiliations}])

def nestDistribution(row):
    files = str2list(row["distribution.url"])
    dates = str2list(row["dateModified"])
    # interleave the two values together
    if(files is None):
        print(f"Missing distribution.contentUrl at row {row.identifier}")
        # raise Exception(f"Missing distribution.contentUrl at row {row}")
    elif(dates is None):
        print(f"Missing distribution.dateModified at row {row.identifier}")
        # raise Exception(f"Missing distribution.dateModified at row {row}")
    elif(len(files) == len(dates)):
        return([{"contentUrl": pair[0], "dateModified": pair[1]} for pair in zip(files, dates)])
    elif(len(dates) == 1):
        return([{"contentUrl": file, "dateModified": dates[0]} for file in files])
    else:
        raise Exception(f"mismatch between files / dateModified at row {row}")

def nestPublisher(database):
    if((database == database) & (database is not None)):
        return({"name": database})

def nestFunding(row):
    ids = str2list(row["funding.identifier"])
    funders = str2list(row["funder.name"])
    # interleave the two values together
    if(ids is None):
        raise Exception(f"Missing funding identifier at row {row}")
    elif(funders is None):
        raise Exception(f"Missing funder name at row {row}")
    elif(len(ids) == len(funders)):
        # paired ids + funders
        return([{"identifier": pair[0], "funder": {"name": pair[1]}} for pair in zip(ids, funders)])
    elif(len(funders) == 1):
        # Assume funder applies to all ids
        return([{"identifier": id, "funder": {"name": funders[0]}} for id in ids])
    else:
        # Assume affiliation applies to all names
        raise Exception(f"mismatch between funding ids / funder names at row {row}")


# Split strings into arrays, as needed
def str2list(value):
    if((value == value) & (value is not None)):
        if((type(value) == int) | (type(value) == float)):
            value = "{:.0f}".format(value)
        return list(map(str.strip, value.split(",")))


    # Pull the publication info
def convertAuthor(authorObj):
    try:
        given = authorObj['given']
    except:
        given = None
    try:
        family = authorObj['family']
    except:
        family = None
    return({'givenName': given,
            'familyName': family
            })

def dateArr2Str(arr):
    if(len(arr) == 3):
        return(f"{str(arr[0])}-{str(arr[1]).zfill(2)}-{str(arr[2]).zfill(2)}")
    if(len(arr) == 2):
        return(f"{str(arr[0])}-{str(arr[1]).zfill(2)}-01")
    if(len(arr) == 1):
        return(f"{str(arr[0])}-01-01")

def getCitation(pmid, isPMID = True):
    if(isPMID):
        ncbi_stub="https://api.ncbi.nlm.nih.gov/lit/ctxp/v1/pubmed/?format=csl&id="
    else:
    # PMCID
        ncbi_stub="https://api.ncbi.nlm.nih.gov/lit/ctxp/v1/pmc/?format=csl&id="
    if(pmid == pmid):
        if(isinstance(pmid, str)):
            pmid_string = pmid
        elif(isinstance(pmid, int) | isinstance(pmid, float)):
            pmid_string = "{:.0f}".format(pmid)
        if(pmid_string):
            res = requests.get(ncbi_stub + pmid_string)
            citation = {}
            if(res.status_code == 200):
                citation_raw = res.json()
                # reformat to schema.org format.
                try:
                    citation['doi'] = citation_raw['DOI']
                except:
                    pass
                citation['@type'] = "ScholarlyArticle"
                citation['pmid'] = citation_raw['PMID']
                try:
                    citation['pmcid'] = citation_raw['PMCID']
                except:
                    pass
                try:
                    citation['doi'] = citation_raw['DOI']
                except:
                    pass
                citation["identifier"] = citation_raw['id']
                citation["issn"] = citation_raw['ISSN']

                citation["author"] = [convertAuthor(
                    author) for author in citation_raw['author']]
                citation["datePublished"] = dateArr2Str(
                    citation_raw['issued']['date-parts'][0])

                try:
                    citation["issueNumber"] = citation_raw['issue']
                except:
                    pass
                citation["journalName"] = citation_raw['container-title']
                citation["journalNameAbbrev"] = citation_raw['container-title-short']
                citation["name"] = citation_raw['title']

                citation["pagination"] = citation_raw['page']
                try:
                    citation["volumeNumber"] = citation_raw['volume']
                except:
                    pass
                citation['url'] = 'https://www.ncbi.nlm.nih.gov/pubmed/?term=' + \
                    citation['pmid']
                return(citation)
            else:
                print(f"WARNING: no citation found for PMID: {pmid_string}")


def getAllCitations(pmids):
    if((pmids == pmids) & (pmids is not None)):
        return([getCitation(id) for id in pmids])


# Script to transform the flat file to one that complies with the schema and includes objects.
def schemaizeMetadata(df, output_name = "NIAID-SysBio-Datasets"):
    # set output file name
    today = date.today().strftime("%Y-%m-%d")
    output = f"{today}_{output_name}.json"

    # define output columns
    output_cols = ["@type", "identifier", "doi", "name", "description", "measurementTechnique", "creator", "sdPublisher", "distribution", "funding", "citation", "license", "species", "infectiousAgent", "healthCondition", "spatialCoverage", "temporalCoverage"]

    # rename any properties, if needed
    df.rename(columns = {"type": "@type", "dataset.doi": "doi"}, inplace=True)

    # Split strings into arrays, as needed
    df["measurementTechnique"] = df.measurementTechnique.apply(str2list)
    df["species"] = df.species.apply(str2list)
    df["infectiousAgent"] = df.infectiousAgent.apply(str2list)
    df["healthCondition"] = df.healthCondition.apply(str2list)
    df["spatialCoverage"] = df.spatialCoverage.apply(str2list)

    # nest related properties into a single object, from the flat columns
    df["sdPublisher"] = df.database.apply(nestPublisher)
    df["creator"] = df.apply(nestAuthor, axis = 1)
    df["distribution"] = df.apply(nestDistribution, axis = 1)
    df["funding"] = df.apply(nestFunding, axis = 1)

    # Pull the publication info
    df["pmid"] = df["publication.pmid"].apply(str2list)
    # df["pmid"] = df["publication.pmid"].apply(str2list)
    df["citation"] = df.pmid.apply(getAllCitations)

    # Check that identifiers are unique
    dupes = df[df.duplicated(["identifier"]) > 0]
    if len(dupes)> 0:
        print(list(dupes.identifier))
        raise Exception(f"{len(dupes)} duplicate identifiers found:")

    df.loc[:, output_cols].to_json(output, orient = "records")

# schemaizeMetadata(df.sample(10))
schemaizeMetadata(df)
