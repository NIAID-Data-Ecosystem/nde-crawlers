from hub.dataload.nde import NDESourceUploader


class DDEUploader(NDESourceUploader):
    name = "dde"
    __metadata__ = {
        "src_meta": {
            "url": "https://discovery.biothings.io/api/dataset/",
            "license_url": "https://creativecommons.org/licenses/by/4.0/",
            "license": "Creative Commons Attribution 4.0 International"
        }
    }