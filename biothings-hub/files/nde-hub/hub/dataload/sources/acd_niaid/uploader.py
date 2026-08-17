from hub.dataload.nde import NDESourceUploader


class AccessClinicalDataUploader(NDESourceUploader):
    name = "acd_niaid"

    __metadata__ = {
        "merger": "merge_struct",
        "merger_kwargs": {"aslistofdict": "", "include": ["includedInDataCatalog"]},
    }
