from hub.dataload.nde import NDESourceUploader


class VEuPathDB_Uploader(NDESourceUploader):
    name = "veupathdb"
    __metadata__ = {
        "merger": "merge_struct",
        "merger_kwargs": {"aslistofdict": "includedInDataCatalog", "include": ["includedInDataCatalog"]},
    }
