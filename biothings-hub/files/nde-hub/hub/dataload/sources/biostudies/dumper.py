import os
import subprocess
from datetime import datetime

import biothings
import biothings.hub.dataload.dumper as dumper
import config

biothings.config_for_app(config)


class Biostudies_Dumper(dumper.WgetDumper):

    SRC_NAME = "biostudies"
    SRC_ROOT_FOLDER = os.path.join(config.DATA_ARCHIVE_ROOT, SRC_NAME)
    # TODO change this into config
    SITEMAP_URLS = config.SITEMAP_URLS

    def set_release(self):
        last_modified_timestamp = os.path.getmtime(self.SITEMAP_URLS)
        # Convert the timestamp to a human-readable format
        self.release = datetime.fromtimestamp(last_modified_timestamp).strftime("%Y-%m-%d")

    def create_todump_list(self, force=False, **kwargs):
        self.set_release()  # so we can generate new_data_folder
        for filename in os.listdir(self.SITEMAP_URLS):
            new_localfile = os.path.join(self.new_data_folder, filename.replace(".txt", ""))
            try:
                current_localfile = os.path.join(self.current_data_folder, filename.replace(".txt", ""))
            except TypeError:
                # current data folder doesn't even exist
                current_localfile = new_localfile

            if force or not os.path.exists(current_localfile):
                remoteurl = os.path.join(self.SITEMAP_URLS, filename)
                self.to_dump.append({"remote": remoteurl, "local": new_localfile})

    def download(self, remoteurl, localfile):
        self.prepare_local_folders(localfile)
        cmdline = "wget -i %s -P %s -N" % (remoteurl, localfile)
        return_code = subprocess.run(["wget", "-i", remoteurl, "-P", localfile, "-N"])
        if return_code.returncode == 0:
            self.logger.info("Success.")
        else:
            self.logger.error("Failed with return code (%s)." % return_code)
