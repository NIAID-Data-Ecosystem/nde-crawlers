import os

import biothings
import biothings.hub.dataload.uploader as uploader
import timeout_decorator
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

    @timeout_decorator.timeout(1800, timeout_exception=TimeoutError)
    def parse_files_with_timeout(self, input_file):
        return list(parse_files(input_file))

    def load_data(self, input_file):
        try:
            return self.parse_files_with_timeout(input_file)
        except TimeoutError:
            self.logger.info("Job timed out, TimeoutError in %s", input_file)
            # return an empty list as BasicStorage expects an iterable
            return []

    @classmethod
    def get_mapping(cls):
        return NDESourceUploader.get_mapping()
