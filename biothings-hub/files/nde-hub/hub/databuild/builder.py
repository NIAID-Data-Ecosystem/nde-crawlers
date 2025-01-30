import biothings.hub.databuild.builder as builder


class NDEDataBuilder(builder.DataBuilder):
    """Merge order for NDE data sources. Highest priority is merged last"""

    def merge_order(self, other_sources):
        self.logger.info("Other sources: %s", other_sources)
        # Priority list of sources to merge from highest to lowest
        priority = ["massive", "ncbi_geo", "lincs", "omicsdi", "veupathdb", "veupath_collections"]
        # Reverse list b/c sources are upserted so highest priority needs to be merged last
        for source in reversed(priority):
            other_sources.append(other_sources.pop(other_sources.index(source)))
        self.logger.info("This is the merge order: %s", other_sources)
        return other_sources
