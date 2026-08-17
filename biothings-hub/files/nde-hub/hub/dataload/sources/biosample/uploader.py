from hub.dataload.nde import NDESourceSampleUploader

from .parser import apply_sex_mapping, load_sex_mapping


class BiosampleUploader(NDESourceSampleUploader):
    name = "biosample"

    _sex_mapping = None

    def post_process(self, doc):
        """Map the sample's reported sex onto our controlled vocabulary."""
        if self._sex_mapping is None:
            self._sex_mapping = load_sex_mapping()
        return apply_sex_mapping(doc, self._sex_mapping)
