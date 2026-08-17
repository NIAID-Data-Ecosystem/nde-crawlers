from hub.dataload.nde import NDESourceUploader


class EmpiarUploader(NDESourceUploader):
    name = "empiar"
    __metadata__ = {
        "merger": "merge_struct",
        "merger_kwargs": {"aslistofdict": "includedInDataCatalog", "include": ["includedInDataCatalog"]},
    }
