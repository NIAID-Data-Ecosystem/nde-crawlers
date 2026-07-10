import re
from collections import defaultdict

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

        projection = {"_id": 1, "doi": 1, "includedInDataCatalog": 1}
        doi_filter = {"doi": {"$exists": True, "$nin": [None, ""]}}

        def prefix_query(source):
            return {**doi_filter, "_id": {"$regex": f"^{re.escape(source)}"}}

        def doi_key(doc):
            doi = doc.get("doi")
            if not doi:
                return None
            if isinstance(doi, list):
                doi = tuple(value for value in doi if value not in (None, ""))
            return doi or None

        def track_doc(docs_by_doi, counts_by_doi, key, doc):
            counts_by_doi[key] += 1
            if len(docs_by_doi[key]) < 2:
                docs_by_doi[key].append(doc)

        def scan_by_prefix(source):
            cursor = collection.find(prefix_query(source), projection, no_cursor_timeout=True).batch_size(1000)
            try:
                yield from cursor
            finally:
                cursor.close()

        self.logger.info(
            f"Scanning documents with duplicate {duplicate} and sources {sources} in {collection_name}"
        )

        source_docs_by_doi = defaultdict(list)
        source_counts_by_doi = defaultdict(int)
        source_doc_count = 0
        for source in sources:
            for doc in scan_by_prefix(source):
                key = doi_key(doc)
                if key is not None:
                    track_doc(source_docs_by_doi, source_counts_by_doi, key, doc)
                    source_doc_count += 1

        self.logger.info(
            f"Found {len(source_docs_by_doi)} DOI groups from {source_doc_count} source records"
        )

        duplicate_docs_by_doi = defaultdict(list)
        duplicate_counts_by_doi = defaultdict(int)
        duplicate_doc_count = 0
        for doc in scan_by_prefix(duplicate):
            key = doi_key(doc)
            if key in source_docs_by_doi:
                track_doc(duplicate_docs_by_doi, duplicate_counts_by_doi, key, doc)
                duplicate_doc_count += 1

        self.logger.info(
            f"Found {len(duplicate_docs_by_doi)} matching DOI groups from {duplicate_doc_count} duplicate records"
        )

        # records to be deleted
        records_to_delete = []
        # bulk operations
        bulk_operations = []
        count = 0

        def flush_writes():
            nonlocal count
            if not bulk_operations:
                return
            collection.bulk_write(bulk_operations)
            count += len(bulk_operations)
            self.logger.info(f"{count} records updated")
            bulk_operations.clear()
            delete_result = collection.delete_many({"_id": {"$in": records_to_delete}})
            self.logger.info(f"Deleted {delete_result.deleted_count} duplicate records in batch")
            records_to_delete.clear()

        group = 0
        skipped_count = 0
        for group, (doi, dupe_docs) in enumerate(duplicate_docs_by_doi.items(), start=1):
            # A single DOI can be shared by the duplicate source (e.g. zenodo), the
            # authoritative source (e.g. dryad/tycho), and occasionally an unrelated
            # source (e.g. figshare re-publishing the dataset under the same DOI).
            # Select exactly one duplicate doc and one source doc to merge, ignoring
            # any other sources, rather than crashing the whole post-merge job.
            source_docs = source_docs_by_doi[doi]
            dupe_count = duplicate_counts_by_doi[doi]
            source_count = source_counts_by_doi[doi]

            if dupe_count != 1 or source_count != 1:
                skipped_count += 1
                sample_ids = [doc["_id"] for doc in [*dupe_docs, *source_docs]]
                self.logger.warning(
                    f"Skipping deduplication group {group}: expected one '{duplicate}' document and one "
                    f"source document, but got {dupe_count} duplicate documents and {source_count} "
                    f"source documents. Sample IDs: {sample_ids}"
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
            records_to_delete.extend(doc["_id"] for doc in dupe_docs)

            if len(bulk_operations) == 1000:
                flush_writes()

        flush_writes()

        self.logger.info(
            f"Merging complete. {group} groups found. Updated {count} records. Skipped {skipped_count} groups"
        )

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

        if "empiar" in source_names:
            self.delete_unmerged_empiar_records()
