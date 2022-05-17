# NIAID Systems Biology Manual importation

This set of code is designed not to be directly imported into the
[NIAID Data Ecosystem metadata API](api.niaid.nih.gov), but rather to prepare
metadata based on a template, to be added to the [Data Discovery Engine NIAID Guide](https://discovery.biothings.io/guide/niaid).

Metadata, in tabular form, are imported and then transformed to an
array of objects to comply with the [NIAID schemas](https://discovery.biothings.io/view/niaid).
The resulting JSON file can be imported into the DDE.
