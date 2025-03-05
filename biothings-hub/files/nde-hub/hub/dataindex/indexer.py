from biothings.hub.dataindex.indexer import Indexer


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
            "filter": ["lowercase", "asciifolding"],
        }

        # Testing standard tokenizer instead of the custom char_tokenizer
        self.es_index_settings["analysis"]["analyzer"]["phrase_suggester"] = {
            "type": "custom",
            "tokenizer": "char_tokenizer",
            "char_filter": ["html_strip"],
            "filter": ["lowercase", "asciifolding", "shingle"],
        }
