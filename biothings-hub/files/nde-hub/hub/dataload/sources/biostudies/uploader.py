import os

import biothings
import biothings.hub.dataload.uploader as uploader
import timeout_decorator
from hub.dataload.nde import NDESourceUploader
from utils.pubtator import standardize_data
from utils.topic_category_helper import add_topic_category
from utils.utils import nde_upload_wrapper

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

    @nde_upload_wrapper
    def load_data(self, input_file):
        try:
            docs = self.parse_files_with_timeout(input_file)
            docs = standardize_data(docs)
            docs = add_topic_category(docs, self.name)
            for doc in docs:
                yield doc
        except TimeoutError:
            self.logger.info("Job timed out, TimeoutError in %s", input_file)
            # return an empty list as BasicStorage expects an iterable
            return []

    @classmethod
    def get_mapping(cls):
        return NDESourceUploader.get_mapping()
