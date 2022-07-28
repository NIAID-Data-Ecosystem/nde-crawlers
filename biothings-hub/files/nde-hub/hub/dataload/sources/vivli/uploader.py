from hub.dataload.nde import NDESourceUploader


class VivliUploader(NDESourceUploader):
    name = "vivli"
    __metadata__ = {
        "src_meta": {
            "url": "https://search.vivli.org/",
            "license_url": "https://vivli.org/resources/vivli-data-use-agreement/"
        }
    }
