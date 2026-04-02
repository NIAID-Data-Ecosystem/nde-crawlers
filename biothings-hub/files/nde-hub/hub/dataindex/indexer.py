from biothings.hub.dataindex.indexer import Indexer
from elasticsearch import AsyncElasticsearch

from .embed import run_embeddings


class NDEIndexer(Indexer):
    def __init__(self, build_doc, indexer_env, index_name):
        super().__init__(build_doc, indexer_env, index_name)

        self.es_index_settings["analysis"]["tokenizer"] = {
            "char_tokenizer": {
                "type": "char_group",
                "tokenize_on_chars": ["whitespace", "\n", ",", "/", "\\", "_", ";", "(", ")"],
            }
        }

        self.es_index_settings["analysis"]["filter"] = {
            "shingle": {"type": "shingle", "min_shingle_size": 2, "max_shingle_size": 3}
        }

        self.es_index_settings["analysis"]["analyzer"]["nde_analyzer"] = {
            "type": "custom",
            "tokenizer": "char_tokenizer",
            "char_filter": ["html_strip"],
            "filter": ["lowercase", "asciifolding", "stemmer"],
        }

        # Testing standard tokenizer instead of the custom char_tokenizer
        self.es_index_settings["analysis"]["analyzer"]["phrase_suggester"] = {
            "type": "custom",
            "tokenizer": "char_tokenizer",
            "char_filter": ["html_strip"],
            "filter": ["lowercase", "asciifolding", "shingle", "stemmer"],
        }

    async def post_index(self, job_manager, *args, **kwargs):
        self.logger.info("Starting post-index embedding for index '%s'", self.es_index_name)

        pinfo = self.pinfo.get_pinfo(
            step="post_index",
            description=f"embedding {self.es_index_name}",
        )

        job = await job_manager.defer_to_process(
            pinfo,
            run_embeddings,
            self.es_client_args.get("hosts"),
            self.es_index_name,
        )
        await job

        self.logger.info("Post-index embedding complete for index '%s'", self.es_index_name)

        self.logger.info("Force merging index '%s' to 1 segment", self.es_index_name)
        async with AsyncElasticsearch(**self.es_client_args) as client:
            await client.indices.forcemerge(index=self.es_index_name, max_num_segments=1, request_timeout=10800)
        self.logger.info("Force merge complete for index '%s'", self.es_index_name)
