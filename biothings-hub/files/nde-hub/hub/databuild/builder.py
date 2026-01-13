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
            "omicsdi",
            "immport",
            "immunespace",
            "veupathdb",
            "veupath_collections",
        ]
        # Reverse list b/c sources are upserted so highest priority needs to be merged last
        for source in reversed(priority):
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
                        "$elemMatch": {"_id": {"$regex": f"^{duplicate}"}},
                        "$elemMatch": {"_id": {"$regex": f"^({'|'.join(sources)})"}},
                    },
                }
            },
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

            if len(result["documents"]) > 2:
                # find all of the duplicate records
                def is_valid_document(doc):
                    if any(doc["_id"].startswith(source) for source in sources):
                        return True
                    if doc["_id"].startswith(duplicate):
                        sdPublisher = doc.get("sdPublisher")
                        if isinstance(sdPublisher, list):
                            return any(
                                any(source in pub.get("name", "").lower() for source in sources) for pub in sdPublisher
                            )
                        elif isinstance(sdPublisher, dict):
                            return any(source in sdPublisher.get("name", "").lower() for source in sources)
                    return False

                result["documents"] = [doc for doc in result["documents"] if is_valid_document(doc)]

            # TODO: If two sources have the same doi this will not work as expected so assert that the length is 2
            # Will also need to handle includedInDataCatalog field if there are more than 2 sources
            assert (
                len(result["documents"]) == 2
            ), f"Expected 2 documents when merging, but got {len(result['documents'])}. Documents: {result['documents']}"

            for doc in result["documents"]:
                if not doc["_id"].startswith(duplicate):
                    _id = doc["_id"]
                    source_includedInDataCatalog = doc.get("includedInDataCatalog")
                else:
                    dupe_includedInDataCatalog = doc.get("includedInDataCatalog")

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

    def identifier_deduplication(self, identifier_source_catalog, source_catalogs):
        """Given a primary source catalog (e.g., Data Discovery Engine) and a list of source catalogs
        (e.g., NCBI GEO), find and match the identifers from the primary source catalog to the _id fields of
        documents from the source catalogs. Merge the resulting documents.

        Args:
            identifier_source_catalog (str): The name of the source identifier catalog to search (e.g., "Data Discovery Engine")
            source_catalogs (list): A list of source catalogs to match against (e.g., ["NCBI GEO", "AccessClinicalData@NIAID"])

        """
        db = mongo.get_target_db()
        collection_name = self.target_backend.target_name
        collection = db[collection_name]

        self.logger.info(
            f"Agregating documents with duplicate {identifier_source_catalog} and sources {source_catalogs} in {collection_name}"
        )

        # Aggregation pipeline to match DDE documents with NCBI GEO documents
        pipeline = [
            # Stage 1: Match documents from Data Discovery Engine
            {"$match": {"includedInDataCatalog.name": identifier_source_catalog}},
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
            # Stage 5: Keep only GEO matches (and drop any non-GEO doc that happened to match)
            {
                "$addFields": {
                    "matching_catalog_docs": {
                        "$filter": {
                            "input": "$matching_docs",
                            "as": "doc",
                            "cond": {
                                "$in": [
                                    "$$doc.includedInDataCatalog.name",
                                    source_catalogs,
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

            donor_doc = result["original_doc"]
            for recipient_doc in result["matching_catalogs"]:
                merged_doc = merge_struct(
                    recipient_doc,
                    donor_doc,
                    aslistofdict="includedInDataCatalog",
                    include=["includedInDataCatalog"],
                )
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

    def post_merge(self, source_names, batch_size, job_manager):
        duplicate = "zenodo"
        sources = ["dryad", "tycho"]
        self.deduplication(sources, duplicate)

        identifier_source_catalog = "Data Discovery Engine"
        source_catalogs = ["NCBI GEO"]
        self.identifier_deduplication(identifier_source_catalog, source_catalogs)
