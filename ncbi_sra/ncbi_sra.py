from pprint import pprint
import requests
import json
import xmltodict


def get_figshare_data(url):

    r = requests.get(url)
    data = xmltodict.parse(r.text)
    json_data = json.dumps(data)
    pprint(json_data)


get_figshare_data(
    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=6678417')
# use LWP::Simple;
# $query = 'chimpanzee[orgn]+AND+biomol+mrna[prop]';

# #assemble the esearch URL
# $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
# $url = $base . "esearch.fcgi?db=nucleotide&term=$query&usehistory=y";

# #post the esearch URL
# $output = get($url);

# #parse WebEnv, QueryKey and Count (# records retrieved)
# $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
# $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
# $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);

# #open output file for writing
# open(OUT, ">chimp.fna") || die "Can't open file!\n";

# #retrieve data in batches of 500
# $retmax = 500;
# for ($retstart = 0; $retstart < $count; $retstart += $retmax) {
#         $efetch_url = $base ."efetch.fcgi?db=nucleotide&WebEnv=$web";
#         $efetch_url .= "&query_key=$key&retstart=$retstart";
#         $efetch_url .= "&retmax=$retmax&rettype=fasta&retmode=text";
#         $efetch_out = get($efetch_url);
#         print OUT "$efetch_out";
# }
# close OUT;

url = ' https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=$query&usehistory=y'
r = requests.get(url)
