import asyncio
import functools

from biothings.hub.dataindex.indexer import Indexer

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

    async def post_index(self, *args, **kwargs):
        es_hosts = self.es_client_args.get("hosts")
        index_name = self.es_index_name
        self.logger.info("Starting post-index embedding for index '%s'", index_name)
        loop = asyncio.get_event_loop()
        await loop.run_in_executor(
            None,
            functools.partial(run_embeddings, es_hosts, index_name, log=self.logger),
        )
        self.logger.info("Post-index embedding complete for index '%s'", index_name)
