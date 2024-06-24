from hub.dataload.nde import NDEFileSystemDumper


class NCBI_Geo_Dumper(NDEFileSystemDumper):
    SRC_NAME = "ncbi_geo"
    SCHEDULE = "0 17 1 1,4,7,10 *"  # 1st day of Jan, Apr, Jul, Oct at 17:00 UTC
