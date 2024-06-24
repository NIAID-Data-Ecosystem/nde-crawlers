from hub.dataload.nde import NDEFileSystemDumper


class OmicsDIDumper(NDEFileSystemDumper):
    SRC_NAME = "omicsdi"
    SCHEDULE = "0 0 1 6,12 *"  # Every 6 months on the 1st day of the month
