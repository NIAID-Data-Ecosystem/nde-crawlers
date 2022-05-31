from datetime import datetime
import logging
import json
import requests


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')


def parse():

    # retrieve entire list of assays with min details
    assay_list_url = f"https://reframedb.org/api/assay_list"
    response = requests.get(assay_list_url)
    assays = json.loads(response.text)

    # save the id of each assay to request later w/ full details
    assay_ids = []
    for assay in assays:
        assay_ids.append(assay['assay_id'])

    # author names too inconsistent, need to manually map each name to affiliation
    author_names = {"Thomas Rogers": "USCD and Scripps Research", "California Institute for Biomedical Research (Calibr), La Jolla, CA": None, 'Brad Spellberg': 'University of Southern California, Los Angeles CA', 'Brian M. Luna': 'University of Southern California, Los Angeles CA', 'Case W. McNamara': 'California Institute for Biomedical Research (Calibr), La Jolla CA', 'FIX:The Global Antibiotic Research and Development Partnership': None, 'Kerstin Gagaring': 'California Institute for Biomedical Research (Calibr), La Jolla CA', 'Malina Bakowski': 'California Institute for Biomedical Research (Calibr), La Jolla CA', 'Michael Green': "University of Massachusetts Medical School", 'Jingxin Wang': 'California Institute for Biomedical Research (Calibr), La Jolla CA', 'Michael Harbut': "California Institute for Biomedical Research (Calibr), La Jolla CA", 'Lindsay Flint': None, 'Yulia Ovechkina': None, 'Alan Chu': None, 'Chen': None, 'Judy Sakanari': 'University of California, San Francisco Dept. of Pharmaceutical Chemistry', 'Melissa S. Love': "California Institute for Biomedical Research (Calibr), La Jolla, CA", 'Martin Redhead': "Head of Quantitative Pharmacology (ExScientia)", 'Joanna Evans': "Center for Discovery and Innovation, Hackensack Meridian Health, Nutley, NJ", 'Prof. Pierre Buffet MD': "INTS, Paris, France", 'Karen Wolff': "California Institute for Biomedical Research (Calibr), La Jolla, CA", 'Changyou Lin': None, 'Cheryl N Miller': "NRC Research Associate, Molecular and Translational Sciences Division, USAMRIID", 'Ryan Choi': "University of Washington, CERID, Seattle, WA", 'Thomas Pietschmann': 'TWINCORE', 'Deli Huang': None, 'Jose L. Lopez-Ribot': 'The University of Texas at San Antonio (UTSA), San Antonio TX', 'Laura Riva': None, 'Malina A. Bakowski': "California Institute for Biomedical Research (Calibr), La Jolla CA", 'Alyssa Manning': None,
                    'Steven De Jonghe': "Rega Institute, KU Leuven, Belgium", 'Robert Bostwick': "(Southern Research, HTS Center), Birmingham, Alabama", 'Monica Kangussu-Marcolino': "Stanford University, StanfordCA", 'Huihui Mou': 'Scripps Research Institute, Jupiter, FL', 'Christopher A. Rice': "University of Georgia", 'Dennis E. Kyle': "University of Georgia", 'Jair Siqueira-Neto': "University of California, San Diego", 'Kristen Johnson': None, 'Emily Chen': None, 'John Malona': "Associate Director (Jnana Therapeutics)", 'Miglianico, Marie': 'TropIQ Health Sciences, The Netherlands', 'Sameera Patel': None, 'Beatrice Cubitt': "Scientific Associate; TSRI, Department of Immunology and Microbiology", 'Erica Penn': "Walter Reed Army Institute of Research, Silver Spring, MD", 'Anjan Debnath': "UCSD, San Diego CA", 'Christopher L. Barratt': "University of Dundee", 'Christina Stallings': None, 'Andrea J. Pruijssers': "Research Assistant Professor of Pediatrics, (Vanderbilt University Medical Center, Nashville, TN, USA)", 'David Nemazee': "Professor, Department of Immunology and Microbiology, The Scripps Research Institute", 'Hui Guo': "California Institute for Biomedical Research (Calibr), La Jolla CA", 'Walter J. Sandoval-Espinola': "Postdoctoral Fellow in Chemistry and Chemical Biology (Harvard University)", 'Sam J Wilson': "MRC University of Glasgow Centre for Virus Research", 'Meredith Davis-Gardner': "The Scripps Research Institute, Jupiter, FL", 'Kenneth Keiler': "Professor of Biochemistry & Molecular Biology, Penn State University", 'Amanda Brown': "Cornell University, Ithaca, NY", 'Tanya Parish': None, 'David G. Russell': "Cornell University, Ithaca, NY", "The Global Antibiotic Research and Development Partnership (GARDP, Geneva, Switzerland": None, 'Franz S. Gruber': 'University of Dundee', 'Yuka Otsuka': "Scripps Research Institute, Jupiter, FL", "Gretchen Ehrenkaufer": "Stanford University, StanfordCA", "Paul Andrews": "University of Dundee", "Feng Wang": "Calibr at Scripps Research", "Chen, Shuibing": None, "Duan, Xiaohua": None, "Dechering, Koen": "TropIQ Health Sciences, The Netherlands", "Linghang Peng": None, "Mitchell V. Hull": "California Institute for Biomedical Research (Calibr), La Jolla CA", "Emily P. Balskus": "Professor of Chemistry and Chemical Biology (Harvard University)", "Jason Wang": None, "Juan C. de la Torre": "TSRI, Department of Immunology and Microbiology", 'Kansas University': None, "Washington University School of Medicine": None, "Yujin Kim": "Postdoctoral Fellow; TSRI, Department of Immunology and Microbiology"}

    count = 0
    for id in assay_ids:
        # ping individual assay url to get full metadata
        individual_assay_url = f"https://reframedb.org/api/assay_details?aid={id}"
        response = requests.get(individual_assay_url)
        # transform the response
        if response.status_code == 200:

            count += 1
            if count % 25 == 0:
                logger.info("Parsed %s records", count)

            metadata = json.loads(response.text)[0]
            output = {"includedInDataCatalog":
                      {"name": "ReframeDB",
                       'versionDate': datetime.today().isoformat(),
                       'url': "https://reframedb.org/"},
                      "@type": "Dataset"
                      }
            if assay_id := metadata.get('assay_id'):
                output['identifier'] = assay_id
                output['_id'] = 'ReframeDB_' + id
            if assay_title := metadata.get('assay_title'):
                output['name'] = assay_title
            if title_short := metadata.get('title_short'):
                output['alternateName'] = title_short
            if authors := metadata.get('authors'):
                authors_list = []
                for key in author_names.keys():
                    if key in authors and author_names[key]:
                        authors_list.append({'name': key,
                                            'affiliation': {'name': author_names[key]}})
                    elif key in authors:
                        authors_list.append({'name': key})
                if len(authors_list) == 0:
                    authors_list.append({'name': authors})
                output['author'] = authors_list
            if summary := metadata.get('summary'):
                output['description'] = summary
            if purpose := metadata.get('purpose'):
                output['description'] += "\nPurpose\n" + purpose
            if protocol := metadata.get('protocol'):
                output['description'] += "\nProtocol\n" + protocol
            if readout := metadata.get('readout'):
                output['description'] += "\nReadout\n" + readout
            if detection_method := metadata.get('detection_method'):
                output['description'] += "\nDetection Method\n" + \
                    detection_method
            if detection_reagents := metadata.get('detection_reagents'):
                output['description'] += "\nDetection Reagents\n" + \
                    detection_reagents
            if components := metadata.get('components'):
                output['description'] += "\nComponents\n" + components
            if drug_conc := metadata.get('drug_conc'):
                output['description'] += "\nDrug Concetration\n" + drug_conc
            if indication := metadata.get('indication'):
                output['healthCondition'] = indication
            if assay_type := metadata.get('assay_type'):
                output['measurementTechnique'] = {'name': assay_type}
            if bibliography := metadata.get('bibliography'):
                output['pmids'] = str(bibliography).split('.')[0]

            yield output

    logger.info("Finished Parsing. Total Records: %s", count)
