from hub.dataload.nde import NDESourceUploader


class GXAUploader(NDESourceUploader):
    name = "gxa"
    __metadata__ = {
        "merger": "merge_struct",
        "merger_kwargs": {"aslistofdict": "includedInDataCatalog", "include": ["includedInDataCatalog"]},
        "src_meta": {
            "url": "https://www.ebi.ac.uk/gxa/",
            "license_url": "https://www.ebi.ac.uk/gxa/licence.html",
        },
    }
