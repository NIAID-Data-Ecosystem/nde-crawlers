import pandas as pd
sheet_id = '1EQbyin1FqdfYmWFglj4lXj6zenugQPAG8SBE2M7osrU'
sheet_name = "Sheet1"
url = f'https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}'
df = pd.read_csv(url)


def get_source_data(source_name):
    col1 = source_name
    col2 = f'nde_{source_name}'
    result = dict(zip(df[col1], df[col2]))
    result = {k: v for k, v in result.items() if not pd.isna(k)
              and not pd.isna(v)}
    return result


{
    "name": "Bioconductor",
    "description": "Bioconductor is a free, open source and open development software project for the analysis and comprehension of genomic data generated by wet lab experiments in molecular biology. It holds a repository of R packages that facilitates rigorous and reproducible analysis of data from current and emerging biological assays. BioConductor delivers releases where a set of packages is published at once and intended for compatibility only with a certain version of R. This is in contrast to CRAN where packages are added continuously with no reference to particular versions of R. Additionally, BioConductor also comes with its own installation tool, BiocManager::install().",
    "schema": get_source_data('bioconductor'),
    "url": "https://bioconductor.org/",
    "identifier": "Bioconductor"
}
