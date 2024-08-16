import os

import biothings.hub.dataload.uploader as uploader
from hub.dataload.nde import NDESourceUploader

from .parser import parse_file_dir


class Biostudies_Uploader(uploader.ParallelizedSourceUploader):
    name = "biostudies"

    def jobs(self):
        jobs = []
        # iterate over each directory in the download directory
        for directory in os.listdir(self.data_folder):
            input_dir = os.path.join(self.data_folder, directory)
            if os.path.isdir(input_dir):
                jobs.append((input_dir,))
        return jobs

    def load_data(self, input_file):
        return parse_file_dir(input_file)

    @classmethod
    def get_mapping(cls):
        return NDESourceUploader.get_mapping()
