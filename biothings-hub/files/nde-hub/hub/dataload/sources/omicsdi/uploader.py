from hub.dataload.nde import NDESourceUploader


class OmicsDIUploader(NDESourceUploader):
    name = "omicsdi"
    __metadata__ = {
        "src_meta": {"url": "https://www.omicsdi.org/search", "license_url": "https://www.ebi.ac.uk/licencing"},
        "merger": "merge_struct",
        "merger_kwargs": {"aslistofdict": "includedInDataCatalog", "include": ["includedInDataCatalog"]},
    }
