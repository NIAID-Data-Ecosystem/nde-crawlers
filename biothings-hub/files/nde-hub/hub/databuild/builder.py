import biothings.hub.databuild.builder as builder
import biothings.utils.mongo as mongo
from biothings.utils.dataload import merge_struct
from pymongo import UpdateOne


class NDEDataBuilder(builder.DataBuilder):
    """Merge order for NDE data sources. Highest priority is merged last"""

    def merge_order(self, other_sources):
        self.logger.info("Other sources: %s", other_sources)
        # Priority list of sources to merge from highest to lowest
        priority = [
            "massive",
            "gse_ncbi_geo",
            "lincs",
            "empiar",
            "omicsdi",
            "immport",
            "immunespace",
            "veupathdb",
            "veupath_collections",
        ]
        # Reverse list b/c sources are upserted so highest priority needs to be merged last
        for source in reversed(priority):
            if source in other_sources:
                other_sources.append(other_sources.pop(other_sources.index(source)))
        self.logger.info("This is the merge order: %s", other_sources)
        return other_sources

    def deduplication(self, sources, duplicate):
        db = mongo.get_target_db()
        collection_name = self.target_backend.target_name
        collection = db[collection_name]

        self.logger.info(f"Agregating documents with duplicate {duplicate} and sources {sources} in {collection_name}")
        pipeline = [
            {"$match": {"doi": {"$ne": None}}},
            {
                "$group": {
                    "_id": "$doi",
                    "documents": {"$push": "$$ROOT"},
                    "count": {"$sum": 1},
                }
            },
            {
                "$match": {
                    "count": {"$gt": 1},
                    "documents": {
                        "$all": [
                            {"$elemMatch": {"_id": {"$regex": f"^{duplicate}"}}},
                            {"$elemMatch": {"_id": {"$regex": f"^({'|'.join(sources)})"}}},
                        ]
                    },
                }
            }
        ]

        # Add allowDiskUse to handle large datasets
        results = collection.aggregate(pipeline, allowDiskUse=True, batchSize=100)

        self.logger.info("Aggregation complete. Merging and deleting duplicate records")

        # records to be deleted
        records_to_delete = []
        # bulk operations
        bulk_operations = []
        count = 0

        # Process results in batches to avoid memory overflow
        batch_count = 0
        for group, result in enumerate(results, start=1):
            batch_count += 1

            # add all records that start with duplicate in _id to records_to_delete
            records_to_delete.extend(
                [doc.get("_id") for doc in result["documents"] if doc["_id"].startswith(duplicate)]
            )

            # A single DOI can be shared by the duplicate source (e.g. zenodo), the
            # authoritative source (e.g. dryad/tycho), and occasionally an unrelated
            # source (e.g. figshare re-publishing the dataset under the same DOI).
            # Select exactly one duplicate doc and one source doc to merge, ignoring
            # any other sources, rather than crashing the whole post-merge job.
            dupe_docs = [doc for doc in result["documents"] if doc["_id"].startswith(duplicate)]
            source_docs = [
                doc for doc in result["documents"] if any(doc["_id"].startswith(source) for source in sources)
            ]

            if len(dupe_docs) != 1 or len(source_docs) != 1:
                self.logger.warning(
                    f"Skipping deduplication group {group}: expected one '{duplicate}' document and one "
                    f"source document, but got {[doc['_id'] for doc in result['documents']]}"
                )
                continue

            _id = source_docs[0]["_id"]
            source_includedInDataCatalog = source_docs[0].get("includedInDataCatalog")
            dupe_includedInDataCatalog = dupe_docs[0].get("includedInDataCatalog")

            # merge the includedInDataCatalog field
            if not isinstance(source_includedInDataCatalog, list):
                source_includedInDataCatalog = [source_includedInDataCatalog]
            if not isinstance(dupe_includedInDataCatalog, list):
                dupe_includedInDataCatalog = [dupe_includedInDataCatalog]

            update_op = UpdateOne(
                {"_id": _id},
                {"$set": {"includedInDataCatalog": source_includedInDataCatalog + dupe_includedInDataCatalog}},
            )
            bulk_operations.append(update_op)

            # Process smaller batches more frequently
            if len(bulk_operations) == 1000:  # Reduced from 10000
                collection.bulk_write(bulk_operations)
                bulk_operations = []
                count += 1000
                self.logger.info(f"{count} records updated")

            # Also process deletion batches to free memory
            if len(records_to_delete) >= 5000:
                delete_result = collection.delete_many({"_id": {"$in": records_to_delete}})
                self.logger.info(f"Deleted {delete_result.deleted_count} duplicate records in batch")
                records_to_delete = []

        if bulk_operations:
            collection.bulk_write(bulk_operations)
            count += len(bulk_operations)
            self.logger.info(f"{count} records updated")

        self.logger.info(
            f"Merging complete. {group} Groups found. Updated {count} records. Deleting remaining duplicate records"
        )

        # Perform final bulk deletion of remaining duplicate records
        if records_to_delete:
            delete_result = collection.delete_many({"_id": {"$in": records_to_delete}})
            self.logger.info(f"Deleted {delete_result.deleted_count} remaining duplicate records")

    def identifier_deduplication(self, identifier_source_catalog, source_catalogs, prefer_matching_catalog_doc=False):
        """Given a source catalog (e.g., Data Discovery Engine) and a list of primary source catalogs
        (e.g., NCBI GEO), find and match the identifers from the source catalog to the _id fields of
        documents from the source catalogs. Merge the resulting documents.

        Args:
            identifier_source_catalog (str): The name of the source identifier catalog to search (e.g., "Data Discovery Engine")
            source_catalogs (list): A list of source catalogs to match against (e.g., ["NCBI GEO", "accessclinicaldata@NIAID"])
            prefer_matching_catalog_doc (bool): When True, keep matching source catalog fields authoritative and
                merge identifier-source catalog provenance into that record.

        """
        db = mongo.get_target_db()
        collection_name = self.target_backend.target_name
        collection = db[collection_name]

        self.logger.info(
            f"Aggregating documents with duplicate {identifier_source_catalog} and sources {source_catalogs} in {collection_name}"
        )

        matching_doc_catalog_names = {
            "$cond": {
                "if": {"$isArray": "$$doc.includedInDataCatalog.name"},
                "then": "$$doc.includedInDataCatalog.name",
                "else": ["$$doc.includedInDataCatalog.name"],
            }
        }

        # Aggregation pipeline to match identifier-source documents with source-catalog documents
        pipeline = [
            # Stage 1: Match unmerged documents from the identifier source catalog
            {
                "$match": {
                    "$and": [
                        {"includedInDataCatalog.name": identifier_source_catalog},
                        {"includedInDataCatalog.name": {"$nin": source_catalogs}},
                    ]
                }
            },
            # Stage 2: keep original doc and ensure identifier is an array
            {
                "$addFields": {
                    "original_doc": "$$ROOT",
                    "identifier_array": {
                        "$cond": {
                            "if": {"$isArray": "$identifier"},
                            "then": "$identifier",
                            "else": ["$identifier"],
                        }
                    },
                }
            },
            # Stage 3: Lowercase identifiers for join
            {
                "$addFields": {
                    "identifier_normalized": {
                        "$map": {
                            "input": "$identifier_array",
                            "as": "id",
                            "in": {"$toLower": {"$toString": "$$id"}},
                        }
                    }
                }
            },
            # Stage 4: Fast lookup by foreignField "_id" (uses the _id index)
            {
                "$lookup": {
                    "from": collection_name,
                    "localField": "identifier_normalized",
                    "foreignField": "_id",
                    "as": "matching_docs",
                }
            },
            # Stage 5: Keep only matches from the requested source catalogs
            {
                "$addFields": {
                    "matching_catalog_docs": {
                        "$filter": {
                            "input": "$matching_docs",
                            "as": "doc",
                            "cond": {
                                "$and": [
                                    {"$ne": ["$$doc._id", "$_id"]},
                                    {
                                        "$gt": [
                                            {
                                                "$size": {
                                                    "$setIntersection": [matching_doc_catalog_names, source_catalogs]
                                                }
                                            },
                                            0,
                                        ]
                                    },
                                ]
                            },
                        }
                    }
                }
            },
            {"$match": {"matching_catalog_docs": {"$ne": []}}},
            # Stage 6: Project useful fields for review
            {
                "$project": {
                    "_id": 1,
                    "original_doc": 1,
                    "matching_count": {"$size": "$matching_catalog_docs"},
                    "matching_ids": "$matching_catalog_docs._id",
                    "matching_catalogs": "$matching_catalog_docs",
                }
            },
        ]

        results = collection.aggregate(pipeline)
        self.logger.info("Aggregation complete. Merging and deleting duplicate records")

        # logger.info(f"Found {len(results)} DDE documents with matches in {source_catalogs}")

        # for doc in results[:2]:
        #     logger.info(f"DDE doc {doc['_id']} matches: {doc['matching_ids']}")
        #     # logger.info(f"original doc: {doc['original_doc']}")
        #     logger.info(f"matching docs: {doc['matching_catalogs']}")

        # records to be deleted
        records_to_delete = []
        # bulk operations
        bulk_operations = []
        count = 0

        for result in results:
            records_to_delete.append(result.get("_id"))

            identifier_doc = result["original_doc"]
            for matching_doc in result["matching_catalogs"]:
                if prefer_matching_catalog_doc:
                    base_doc = matching_doc
                    provenance_doc = identifier_doc
                else:
                    base_doc = identifier_doc
                    provenance_doc = matching_doc
                merged_doc = merge_struct(
                    base_doc,
                    provenance_doc,
                    aslistofdict="includedInDataCatalog",
                    include=["includedInDataCatalog"],
                )
                merged_doc["_id"] = matching_doc["_id"]  # keep the matching source catalog id
                bulk_operations.append(
                    UpdateOne(
                        {"_id": merged_doc["_id"]},
                        {"$set": merged_doc},
                        upsert=True,
                    )
                )
                count += 1

                if len(bulk_operations) == 1000:
                    collection.bulk_write(bulk_operations)
                    bulk_operations = []
                    count += 1000
                    self.logger.info(f"{count} records updated")

        if bulk_operations:
            collection.bulk_write(bulk_operations)
            self.logger.info(f"Final batch: {len(bulk_operations)} records updated")

        if records_to_delete:
            delete_result = collection.delete_many({"_id": {"$in": records_to_delete}})
            self.logger.info(f"Deleted {delete_result.deleted_count} duplicate records")

    def delete_unmerged_empiar_records(self):
        """Remove direct EMPIAR records that did not merge with OMICS-DI."""
        db = mongo.get_target_db()
        collection_name = self.target_backend.target_name
        collection = db[collection_name]

        delete_result = collection.delete_many(
            {
                "includedInDataCatalog.name": "Electron Microscopy Public Image Archive",
                "$nor": [{"includedInDataCatalog.name": "Omics Discovery Index (OmicsDI)"}],
            }
        )
        self.logger.info(
            "Deleted %d direct-only EMPIAR records that did not merge with OMICS-DI",
            delete_result.deleted_count,
        )

    def post_merge(self, source_names, batch_size, job_manager):
        duplicate = "zenodo"
        sources = ["dryad", "tycho"]
        self.deduplication(sources, duplicate)

        identifier_source_catalog = "Data Discovery Engine"
        source_catalogs = ["NCBI GEO", "MassIVE", "NCBI BioProject", "Protein Data Bank"]
        self.identifier_deduplication(identifier_source_catalog, source_catalogs)

        self.identifier_deduplication(
            "ProteomeXchange",
            ["MassIVE"],
            prefer_matching_catalog_doc=True,
        )

        # CEIRR reagents can also be cataloged in BEI Resources. The CEIRR
        # crawler emits BEI identifiers in the matching BEI _id form so this
        # can merge CEIRR catalog provenance into the BEI record and remove
        # the duplicate CEIRR document.
        self.identifier_deduplication(
            "Centers of Excellence for Influenza Research and Response (CEIRR) Resources",
            ["BEI Resources"],
            prefer_matching_catalog_doc=True,
        )

        if "empiar" in source_names:
            self.delete_unmerged_empiar_records()
