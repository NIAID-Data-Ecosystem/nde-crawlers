import re
import pandas as pd
import datetime
import logging
import requests


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('nde-logger')


def parse():
    logger.info('Retrieving package metadata from Bioconductor')
    url = "https://bioconductor.org/packages/release/bioc/VIEWS"
    data = requests.get(url).text
    pat = re.compile(
        r"^([^\s][^:]*): (.+?)\s*(?=^[^\s][^:]*:|\Z)", flags=re.S | re.M
    )
    package_dict_list = []
    for chunk in data.split("\n\n"):
        if chunk:
            package_dict_list.append(dict(pat.findall(chunk)))

    logger.info('Retrieving Download Statistics from Bioconductor')
    download_stats_df = pd.read_csv(
        'http://bioconductor.org/packages/stats/bioc/bioc_pkg_stats.tab', sep='\t')
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

    count = 0
    for metadata in package_dict_list:
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
            output['identifier'] = identifier
            output['url'] = 'https://bioconductor.org/packages/release/bioc/html/' + \
                identifier+'.html'
            output['doi'] = f'10.18129/B9.bioc.{identifier}'

            downloads = [d['Nb_of_downloads']
                         for d in download_stats_dicts if d['Package'] == identifier][0]
            output['aggregateRating'] = {
                'ratingValue': downloads, 'reviewAspect': 'Downloads in the last 12 months'}

        if version := metadata.get('Version'):
            output['softwareVersion'] = version

        if depends := metadata.get('Depends'):
            depends = [d.replace('\n', '').strip().replace('    ', '')
                       for d in depends.split(',') if d]
            output['softwareRequirements'] = depends

        if license := metadata.get('License'):
            output['license'] = license
        if title := metadata.get('Title'):
            title = title.replace('\n', '').strip().replace('    ', '')
            output['name'] = title
        if description := metadata.get('Description'):
            description = description.replace(
                '\n', '').strip().replace('    ', '')
            output['description'] = description
        if bioc_views := metadata.get('biocViews'):
            bioc_views = [d.replace('\n', '').strip().replace('    ', '')
                          for d in bioc_views.split(',') if d]
            output['keywords'] = bioc_views

        author_list = []
        if authors := metadata.get('Author'):
            authors = re.sub("[\(\[].*?[\)\]]", "", authors)
            authors = [d.replace('\n', '').strip().replace('    ', '')
                       for d in authors.split(',') if d]
            for author in authors:
                if ' and ' in author:
                    url = None
                    one = author.split(' and ')[0]
                    if one.find('<') != -1:
                        url = one[one.find('<')+1:one.find('>')]
                        one = re.sub("[<].*?[>]", "", one)
                        if url:
                            author_list.append({'name': one, 'url': url})
                        else:
                            author_list.append({'name': one})
                    url = None
                    two = author.split(' and ')[1]
                    if two.find('<') != -1:
                        url = two[two.find('<')+1:two.find('>')]
                        two = re.sub("[<].*?[>]", "", two)
                    if url:
                        author_list.append({'name': two, 'url': url})
                    else:
                        author_list.append({'name': two})
                    url = None
                else:
                    if author.find('<') != -1:
                        url = author[author.find('<')+1:author.find('>')]
                    if url:
                        author_list.append({'name': author, 'url': url})
                    else:
                        author_list.append({'name': author})

        if maintainers := metadata.get('Maintainer'):
            maintainers = re.sub("[\(\[].*?[\)\]]", "", maintainers)
            maintainers = [d.replace('\n', '').strip().replace('    ', '')
                           for d in maintainers.split(',') if d]
            for author in maintainers:
                if ' and ' in author:
                    url = None
                    one = author.split(' and ')[0]
                    if one.find('<') != -1:
                        url = one[one.find('<')+1:one.find('>')]
                        one = re.sub("[<].*?[>]", "", one)
                    if url and 'hpages.on.github' not in url:
                        author_list.append({'name': one, 'url': url})
                    else:
                        author_list.append({'name': one})
                    url = None
                    two = author.split(' and ')[1]
                    if two.find('<') != -1:
                        url = two[two.find('<')+1:two.find('>')]
                        two = re.sub("[<].*?[>]", "", two)
                    if url and 'hpages.on.github' not in url:
                        author_list.append({'name': two, 'url': url})
                    else:
                        author_list.append({'name': two})
                    url = None
                else:
                    if author.find('<') != -1:
                        url = author[author.find('<')+1:author.find('>')]
                        author = re.sub("[<].*?[>]", "", author)
                    if url and 'hpages.on.github' not in url:
                        author_list.append({'name': author, 'url': url})
                    else:
                        author_list.append({'name': author})
        for author_dict in author_list:
            if author_dict['name'] in [a['name'] for a in author_list]:
                author_list.remove(author_dict)
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

        software_help = []
        if vignettes := metadata.get('vignettes'):
            vignettes = [d.replace('\n', '').strip().replace('    ', '')
                         for d in vignettes.split(',') if d]
            if vignette_titles := metadata.get('vignetteTitles'):
                vignette_titles = [d.replace('\n', '').strip().replace('    ', '')
                                   for d in vignette_titles.split(',') if d]
                for vignette, vignette_title in zip(vignettes, vignette_titles):
                    software_help.append({
                        'name': vignette_title,
                        'url': 'https://www.bioconductor.org/packages/release/bioc/'+vignette})
        if len(software_help):
            output['softwareHelp'] = software_help

        based_on = []
        if imports := metadata.get('Imports'):
            imports = [d.replace('\n', '').strip().replace('    ', '')
                       for d in imports.split(',') if d]
            for pkg_import in imports:
                pkg_import = re.sub("[()].*?[)]", "", pkg_import)
                pkg_import = pkg_import.strip()
                based_on.append(
                    {
                        'name': pkg_import,
                        'identifier': f'Bioconductor_{pkg_import}',
                        'url': f'https://www.bioconductor.org/packages/release/bioc/html/{pkg_import}.html'
                    }
                )
        if len(based_on):
            output['isBasedOn'] = based_on

        is_enhanced_by = []
        if enhances := metadata.get('Enhances'):
            enhances = [d.replace('\n', '').strip().replace('    ', '')
                        for d in enhances.split(',') if d]
            for pkg_enhance in enhances:
                pkg_enhance = re.sub("[()].*?[)]", "", pkg_enhance)
                pkg_enhance = pkg_enhance.strip()
                is_enhanced_by.append({'identifier': pkg_enhance})
        if len(is_enhanced_by):
            output['softwareAddOn'] = is_enhanced_by

        is_dependent_on = []
        if depends_on_me := metadata.get('dependsOnMe'):
            depends_on_me = [d.replace('\n', '').strip().replace('    ', '')
                             for d in depends_on_me.split(',') if d]
            for pkg_dependency in depends_on_me:
                pkg_dependency = re.sub("[()].*?[)]", "", pkg_dependency)
                pkg_dependency = pkg_dependency.strip()
                is_dependent_on.append({'identifier': pkg_dependency})
        if len(is_dependent_on):
            output['isBasisFor'] = is_dependent_on

        is_related_to = []
        if suggests := metadata.get('Suggests'):
            suggests = [d.replace('\n', '').strip().replace('    ', '')
                        for d in suggests.split(',') if d]
            for pkg_suggest in suggests:
                pkg_suggest = re.sub("[()].*?[)]", "", pkg_suggest)
                pkg_suggest = pkg_suggest.strip()
                is_related_to.append(
                    {
                        'name': pkg_suggest,
                        'identifier': f'Bioconductor_{pkg_suggest}',
                        'url': f'https://www.bioconductor.org/packages/release/bioc/html/{pkg_suggest}.html'
                    }
                )
        if depends := metadata.get('dependsOnMe'):
            depends = [d.replace('\n', '').strip().replace('    ', '')
                       for d in depends.split(',') if d]
            for pkg_dependency in depends:
                pkg_dependency = re.sub("[()].*?[)]", "", pkg_dependency)
                pkg_dependency = pkg_dependency.strip()
                is_related_to.append(
                    {
                        'name': pkg_dependency,
                        'identifier': f'Bioconductor_{pkg_dependency}',
                        'url': f'https://www.bioconductor.org/packages/release/bioc/html/{pkg_dependency}.html'
                    }
                )
        if suggests_me := metadata.get('suggestsMe'):
            suggests_me = [d.replace('\n', '').strip().replace('    ', '')
                           for d in suggests_me.split(',') if d]
            for pkg_suggest in suggests_me:
                pkg_suggest = re.sub("[()].*?[)]", "", pkg_suggest)
                pkg_suggest = pkg_suggest.strip()
                is_related_to.append(
                    {
                        'name': pkg_suggest,
                        'identifier': f'Bioconductor_{pkg_suggest}',
                        'url': f'https://www.bioconductor.org/packages/release/bioc/html/{pkg_suggest}.html'
                    }
                )
        if linking_to := metadata.get('LinkingTo'):
            linking_to = [d.replace('\n', '').strip().replace('    ', '')
                          for d in linking_to.split(',') if d]
            for pkg_link in linking_to:
                pkg_link = re.sub("[()].*?[)]", "", pkg_link)
                pkg_link = pkg_link.strip()
                is_related_to.append(
                    {
                        'name': pkg_link,
                        'identifier': f'Bioconductor_{pkg_link}',
                        'url': f'https://www.bioconductor.org/packages/release/bioc/html/{pkg_link}.html'
                    }
                )
        if len(is_related_to):
            output['isRelatedTo'] = is_related_to

        git_urls = []
        if github_url := metadata.get('URL'):
            github_url = [d.replace('\n', '').strip().replace('    ', '')
                          for d in github_url.split(',') if d]
            for url in github_url:
                url = re.sub("[(].*?[)]", "", url)
                if 'https://github.com/' in url and not 'https://github.com/compbiomed/TBSignatureProfiler' in url:
                    git_urls.append(url)
        if len(git_urls):
            output['codeRepository'] = git_urls

        if bug_reports := metadata.get('BugReports'):
            output['discussionUrl'] = bug_reports

        if archs := metadata.get('Archs'):
            output['processorRequirements'] = archs

        yield output

    logger.info("Finished Parsing. Total Records: %s", count)
