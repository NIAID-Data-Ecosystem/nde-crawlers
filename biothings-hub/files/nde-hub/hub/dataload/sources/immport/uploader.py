from hub.dataload.nde import NDESourceUploader


class ImmPortUploader(NDESourceUploader):
    name = "immport"
    __metadata__ = {
        "src_meta": {
            "url": "https://www.immport.org/shared/home",
            "license_url": "https://docs.immport.org/home/agreement/"
        }
    }

