# this file is more or less just copied over from other projects
import logging
import os
import os.path
from functools import partial

import biothings
from biothings.hub import HubServer
from biothings.utils.version import set_versions

# shut some mouths
logging.getLogger("botocore").setLevel(logging.ERROR)
logging.getLogger("boto3").setLevel(logging.ERROR)
logging.getLogger("s3transfer").setLevel(logging.ERROR)
logging.getLogger("urllib3").setLevel(logging.ERROR)


app_folder, _src = os.path.split(os.path.split(os.path.abspath(__file__))[0])
import config

set_versions(config, app_folder)
biothings.config_for_app(config)
logging = config.logger
import hub.dataload.sources


class NDEHubServer(HubServer):
    def configure_build_manager(self):
        import biothings.hub.databuild.builder as builder
        from hub.databuild.builder import NDEDataBuilder

        # set specific managers
        build_manager = builder.BuilderManager(builder_class=NDEDataBuilder, job_manager=self.managers["job_manager"])
        build_manager.configure()
        self.managers["build_manager"] = build_manager
        self.logger.info("Using custom builder %s" % NDEDataBuilder)


server = NDEHubServer(hub.dataload.sources, name="BioThings Studio for NDE")

# import ptvsd

# Allow other computers to attach to ptvsd at this IP address and port.
# ptvsd.enable_attach(address=('1.2.3.4', 3000), redirect_output=True)

# Pause the program until a remote debugger is attached
# ptvsd.wait_for_attach()

if __name__ == "__main__":
    # vanilla or as a launcher of an API
    import glob
    import re
    import sys
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-a", "--api-folder", help="API folder to run", dest="api_folder")
    (options, args) = parser.parse_args()
    if options.api_folder:
        api_folder = os.path.abspath(options.api_folder)
        logging.info("Lauching server from API located in: %s", api_folder)
        origwd = os.path.abspath(os.path.curdir)
        # assuming a bin/hub.py module in this folder
        os.chdir(api_folder)
        assert "bin" in os.listdir(), "Can't find 'bin' folder containing hub.py"
        scripts = glob.glob(os.path.join("bin", "*.py"))
        print(scripts)
        startup = None
        if len(scripts) == 1:
            startup = scripts.pop()
        else:
            if "bin/hub.py" in scripts:
                startup = scripts[scripts.index("bin/hub.py")]
            else:
                logging.error(
                    "Found more than one startup scripts, none of them named hub.py, " "don't know which to choose: %s",
                    scripts,
                )
                sys.exit(1)
        logging.info("Found startup script '%s'", startup)
        strmod = re.sub(".py$", "", startup).replace("/", ".")
        import importlib

        mod = importlib.import_module(strmod)
        from biothings.hub import HubServer

        # try to locate a hub server instance
        for name in dir(mod):
            server = getattr(mod, name)
            if issubclass(server.__class__, HubServer):
                logging.info("Found hub server: %s", server)
                # replace sources, dynamic discovery
                if os.path.exists("hub/dataload/sources"):
                    server.source_list = "hub/dataload/sources"
                    logging.info("Auto-discovering sources in 'hub/dataload/sources'")
                else:
                    logging.warning(
                        "Couldn't locate sources folder (expecting 'hub/dataload/sources'), keep those defined by the API"
                    )
                break
    else:
        logging.info("Runing vanilla studio")

    logging.info("Hub DB backend: %s", biothings.config.HUB_DB_BACKEND)
    logging.info("Hub database: %s", biothings.config.DATA_HUB_DB_DATABASE)

    server.start()
