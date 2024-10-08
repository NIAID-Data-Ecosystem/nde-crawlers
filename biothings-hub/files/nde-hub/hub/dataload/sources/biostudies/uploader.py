import os

import biothings
import biothings.hub.dataload.uploader as uploader
from hub.dataload.nde import NDESourceUploader

from .parser import parse_files


class Biostudies_Uploader(uploader.ParallelizedSourceUploader):
    name = "biostudies"
    storage_class = biothings.utils.storage.IgnoreDuplicatedStorage
    MAX_PARALLEL_UPLOAD = 3

    def jobs(self):
        jobs = []
        # iterate over each directory in the download directory
        for file in os.listdir(self.data_folder):
            accno_file = os.path.join(self.data_folder, file)
            jobs.append((accno_file,))
        return jobs

    def load_data(self, input_file):
        return parse_files(input_file)

    @classmethod
    def get_mapping(cls):
        return NDESourceUploader.get_mapping()
