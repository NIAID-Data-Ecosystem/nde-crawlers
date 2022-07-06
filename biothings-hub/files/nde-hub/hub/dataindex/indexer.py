from biothings.hub.dataindex.indexer import Indexer

class NDEIndexer(Indexer):
    def __init__(self, build_doc, indexer_env, index_name):
        super().__init__(build_doc, indexer_env, index_name)

        self.es_index_settings["analysis"]["tokenizer"] = {
            "char_tokenizer": {
                "type": "char_group",
                "tokenize_on_chars": [
                    "whitespace",
                    "\n",
                    ",",
                    "/",
                    "\\",
                    "_",
                    ";"
                ]
            }
        }
        
        self.es_index_settings["analysis"]["analyzer"]["nde_analyzer"] = {
            "type": "custom",
            "tokenizer": "char_tokenizer",
            "char_filter": [
                "html_strip"
            ],
            "filter": [
                "lowercase",
                "asciifolding"
            ]
        }







