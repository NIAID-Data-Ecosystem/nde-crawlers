# NIAID Systems Biology Manual importation

This set of code is designed not to be directly imported into the
[NIAID Data Ecosystem metadata API](api.niaid.nih.gov), but rather to prepare
metadata based on a template, to be added to the [Data Discovery Engine NIAID Guide](https://discovery.biothings.io/guide/niaid).

Metadata, in tabular form, are imported and then transformed to an
array of objects to comply with the [NIAID schemas](https://discovery.biothings.io/view/niaid).
The resulting JSON file can be imported into the DDE.

Steps to execute:
1. Update the [Google sheets](https://docs.google.com/spreadsheets/d/12fm2Qkh39Ben7QS5QtIElVlaYwYY8RUExaFX9EX6iAg/edit#gid=91274120) with metadata.
2. If needed, update the [sheet id](https://github.com/NIAID-Data-Ecosystem/nde-crawlers/blob/bf299e0a12c4f51dc4c5dc549bcf88df3923fd90/niaid_sysbio_manual/clean_sysbio_datasets.py#L8) in the cleanup code
3. Run [clean_sysbio_datasets.py](https://github.com/NIAID-Data-Ecosystem/nde-crawlers/blob/main/niaid_sysbio_manual/clean_sysbio_datasets.py)
4. Upload each of the .json outputs to the [DDE SysBio Dataset bulk uploader](https://discovery.biothings.io/guide/niaid)

Assuming the upload is successful, the results will immediately be viewable in the [DDE API](https://discovery.biothings.io/dataset?template=niaid:dataset). To import into the DDE, you will need to run the [NDE/DDE crawler]](https://github.com/NIAID-Data-Ecosystem/nde-crawlers/tree/main/dde_niaid) and publish a new release.
