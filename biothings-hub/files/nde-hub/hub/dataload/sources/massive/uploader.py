from hub.dataload.nde import NDESourceUploader


class MassiveUploader(NDESourceUploader):
    name = "massive"

    __metadata__ = {
        "merger": "merge_struct",
        "merger_kwargs": {"aslistofdict": "includedInDataCatalog", "include": ["includedInDataCatalog"]},
    }
