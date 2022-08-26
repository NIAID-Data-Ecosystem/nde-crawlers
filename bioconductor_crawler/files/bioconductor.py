import pandas as pd
import rpy2
import datetime
import logging

from rpy2.robjects.packages import importr
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')

rpy2_logger.setLevel(logging.ERROR)   # will display errors, but not warnings


def parse():

    # Install the required R packages
    logger.info('Installing the required R packages')
    rpy2.robjects.r('''
        if (!require("BiocManager", quietly = TRUE))
            install.packages("BiocManager")

        BiocManager::install("BiocPkgTools", ask = FALSE, quietly = TRUE)
        ''')

    # Import the top level R package
    logger.info('Importing the top level R package')
    importr('BiocManager')

    # Import the r package that contains a table of metadata
    logger.info('Importing the R package that contains a table of metadata')
    biocPkgTools = importr('BiocPkgTools')

    # Package metadata table
    biocPkgList = biocPkgTools.biocPkgList()

    # Download statistics table
    download_stats = biocPkgTools.biocDownloadStats()

    # Convert the download statistics table to a pandas dataframe, aggregate the data, and convert it to a list of dictionaries
    logger.info('Retrieving download statistics')
    download_stats_df = pd.DataFrame(download_stats)
    download_stats_columns = list(download_stats.colnames)
    download_stats_df = download_stats_df.T
    download_stats_df.columns = download_stats_columns
    download_stats_df = download_stats_df[download_stats_df['Month'] != 'all']
    download_stats_df['date'] = download_stats_df['Month'].astype(
        str) + '/' + download_stats_df['Year'].astype(str)
    download_stats_df['date'] = pd.to_datetime(
        download_stats_df['date'], format='%b/%Y')
    download_stats_df_last_12_months = download_stats_df[download_stats_df['date']
                                                         >= datetime.datetime.now() - datetime.timedelta(days=365)]
    download_stats_df_last_12_months = download_stats_df_last_12_months[
        download_stats_df_last_12_months['date'] <= datetime.datetime.now()]
    aggregation_functions = {'Nb_of_downloads': 'sum'}
    download_stats_df_last_12_months = download_stats_df_last_12_months.groupby(
        download_stats_df_last_12_months['Package']).aggregate(aggregation_functions)
    download_stats_df_last_12_months = download_stats_df_last_12_months.reset_index()
    download_stats_dicts = download_stats_df_last_12_months.to_dict(
        orient='records')

    # Convert the package metadata table to a pandas dataframe, convert it to a list of dictionaries
    logger.info('Retrieving package metadata')
    biocPkgList_df = pd.DataFrame(biocPkgList)
    columns = list(biocPkgList.colnames)
    biocPkgList_df = biocPkgList_df.T
    biocPkgList_df.columns = columns
    biocPkgList_df = biocPkgList_df.applymap(lambda x: None if isinstance(
        x, rpy2.rinterface_lib.sexp.NACharacterType) or isinstance(
        x[0], rpy2.rinterface_lib.sexp.NACharacterType) else x)
    biocPkgList_df = biocPkgList_df.applymap(lambda x: list(x) if isinstance(
        x, rpy2.robjects.vectors.StrVector) else x)
    biocPkgList_dict = biocPkgList_df.to_dict('records')

    count = 0
    for metadata in biocPkgList_dict:
        count += 1
        if count % 100 == 0:
            logger.info(f'Retrieved {count} package metadata')

        output = {'includedInDataCatalog': {
            '@type': 'ComputationalTool',
            'name': 'Bioconductor',
            'url': 'https://bioconductor.org/',
            'versionDate': datetime.date.today().isoformat()
        },
            '@type': 'ComputationalTool',
            'programmingLanguage': 'R',
            'sdPublisher': {'name': 'Bioconductor'},
            'applicationSuite': 'Bioconductor',
        }

        if identifier := metadata.get('Package'):
            output['_id'] = 'Bioconductor_' + identifier
            output['identifer'] = identifier
            output['url'] = 'https://bioconductor.org/packages/release/bioc/html/' + \
                identifier+'.html'

            downloads = [d['Nb_of_downloads']
                         for d in download_stats_dicts if d['Package'] == identifier][0]
            output['aggregateRating'] = {
                'ratingValue': downloads, 'reviewAspect': 'Downloads in the last 12 months'}

        if version := metadata.get('Version'):
            output['softwareVersion'] = version

        if depends := metadata.get('Depends'):
            output['softwareRequirements'] = depends

        if license := metadata.get('License'):
            output['license'] = license
        if title := metadata.get('Title'):
            output['name'] = title
        if description := metadata.get('Description'):
            output['description'] = description
        if bioc_views := metadata.get('biocViews'):
            output['keywords'] = bioc_views

        author_list = []
        if authors := metadata.get('Author'):
            if isinstance(authors, list):
                for author in authors:
                    author_list.append({'name': author})
            else:
                author_list.append({'name': authors})
        if maintainers := metadata.get('Maintainer'):
            if isinstance(maintainers, list):
                for maintainer in maintainers:
                    maintainer_name = maintainer.split(' ')[:2]
                    author_dict = {'name': ' '.join(maintainer_name)}
                    if author_dict not in author_list:
                        author_list.append(author_dict)
            else:
                maintainer_name = maintainers.split(' ')[:2]
                author_dict = {'name': ' '.join(maintainer_name)}
                if author_dict not in author_list:
                    author_list.append(author_dict)
        if len(author_list):
            output['author'] = author_list

        if date_modified := metadata.get('git_last_commit_date'):
            output['dateModified'] = date_modified
        if date_published := metadata.get('Date/Publication'):
            output['datePublished'] = date_published

        download_urls = []
        if download_url := metadata.get('source.ver'):
            download_urls.append({
                'name': 'https://www.bioconductor.org/packages/release/bioc/'+download_url})
        if download_url := metadata.get('win.binary.ver'):
            download_urls.append({
                'name': 'https://www.bioconductor.org/packages/release/bioc/'+download_url})
        if download_url := metadata.get('mac.binary.ver'):
            download_urls.append({
                'name': 'https://www.bioconductor.org/packages/release/bioc/'+download_url})
        if len(download_urls):
            output['downloadUrl'] = download_urls

        if vignettes := metadata.get('vignettes'):
            if vignette_titles := metadata.get('vignetteTitles'):
                if isinstance(vignette_titles, list) and isinstance(vignettes, list):
                    for vignette, vignette_title in zip(vignettes, vignette_titles):
                        output['softwareHelp'] = {
                            'name': vignette_title,
                            'url': 'https://www.bioconductor.org/packages/release/bioc/'+vignette}
                else:
                    output['softwareHelp'] = {
                        'name': vignette_titles,
                        'url': 'https://www.bioconductor.org/packages/release/bioc/'+vignettes}

        if imports := metadata.get('Imports'):
            based_on = []
            if isinstance(imports, list):
                for pkg_import in imports:
                    based_on.append({'identifier': pkg_import})
            else:
                based_on.append({'identifier': imports})
            if len(based_on):
                output['isBasedOn'] = based_on

        if enhances := metadata.get('Enhances'):
            is_enhanced_by = []
            if isinstance(enhances, list):
                for pkg_enhance in enhances:
                    is_enhanced_by.append({'identifier': pkg_enhance})
            else:
                is_enhanced_by.append({'identifier': enhances})
            if len(is_enhanced_by):
                output['softwareAddOn'] = is_enhanced_by

        if depends_on_me := metadata.get('dependsOnMe'):
            is_dependent_on = []
            if isinstance(depends_on_me, list):
                for pkg_dependency in depends_on_me:
                    is_dependent_on.append({'identifier': pkg_dependency})
            else:
                is_dependent_on.append({'identifier': depends_on_me})
            if len(is_dependent_on):
                output['isBasisFor'] = is_dependent_on

        is_related_to = []
        if suggests := metadata.get('Suggests'):
            if isinstance(suggests, list):
                for pkg_suggest in suggests:
                    is_related_to.append({'identifier': pkg_suggest})
            else:
                is_related_to.append({'identifier': suggests})
        if depends := metadata.get('dependsOnMe'):
            if isinstance(depends, list):
                for pkg_dependency in depends:
                    is_related_to.append({'identifier': pkg_dependency})
            else:
                is_related_to.append({'identifier': depends})
        if suggests_me := metadata.get('suggestsMe'):
            if isinstance(suggests_me, list):
                for pkg_suggest in suggests_me:
                    is_related_to.append({'identifier': pkg_suggest})
            else:
                is_related_to.append({'identifier': suggests_me})
        if linking_to := metadata.get('LinkingTo'):
            if isinstance(linking_to, list):
                for pkg_link in linking_to:
                    is_related_to.append({'identifier': pkg_link})
            else:
                is_related_to.append({'identifier': linking_to})
        if len(is_related_to):
            output['isRelatedTo'] = is_related_to

        if github_url := metadata.get('URL'):
            if 'github' in github_url:
                output['codeRepository'] = github_url

        if bug_reports := metadata.get('BugReports'):
            output['discussionUrl'] = bug_reports

        if archs := metadata.get('Archs'):
            output['processorRequirements'] = archs

        yield output

    logger.info("Finished Parsing. Total Records: %s", count)
