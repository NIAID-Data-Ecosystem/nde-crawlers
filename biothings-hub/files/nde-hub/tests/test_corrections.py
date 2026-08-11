import itertools
import importlib.util
import logging
import sys
import tempfile
import types
import unittest
from itertools import islice
from pathlib import Path
from unittest.mock import patch


if not hasattr(itertools, "batched"):

    def batched(iterable, n):
        iterator = iter(iterable)
        while batch := tuple(islice(iterator, n)):
            yield batch

    itertools.batched = batched


HUB_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = HUB_ROOT / "utils" / "corrections.py"

config = sys.modules.get("config")
if config is None:
    config = types.ModuleType("config")
    sys.modules["config"] = config
config.logger = logging.getLogger("test_corrections")
config.token = ""

package_name = "corrections_test_utils"
package = types.ModuleType(package_name)
package.__path__ = [str(HUB_ROOT / "utils")]
sys.modules[package_name] = package

common = types.ModuleType(f"{package_name}.common")
common.as_list = lambda value: [] if value is None else value if isinstance(value, list) else [value]
sys.modules[common.__name__] = common

spec = importlib.util.spec_from_file_location(f"{package_name}.corrections", MODULE_PATH)
corrections = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = corrections
spec.loader.exec_module(corrections)


SAMPLE_INDEX = {
    "by_id": {
        "record": [
            {
                "name": "example",
                "organizations": [{"name": "Example Organization"}],
                "approved": True,
            }
        ]
    },
    "by_funding": [],
}


class CorrectionsCacheTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.cache_folder_patch = patch.object(
            corrections.config,
            "CACHE_FOLDER",
            self.temp_dir.name,
            create=True,
        )
        self.cache_folder_patch.start()
        corrections._corrections_cache = None

    def tearDown(self):
        corrections._corrections_cache = None
        self.cache_folder_patch.stop()
        self.temp_dir.cleanup()

    def test_separate_process_cache_loads_share_one_github_build(self):
        with patch.object(corrections, "_build_corrections_index", return_value=SAMPLE_INDEX) as build:
            first = corrections.get_corrections_index()
            corrections._corrections_cache = None
            second = corrections.get_corrections_index()

        self.assertEqual(first, SAMPLE_INDEX)
        self.assertEqual(second, SAMPLE_INDEX)
        build.assert_called_once_with()

    def test_failed_build_is_not_retried_for_every_document(self):
        with patch.object(corrections, "_load_or_build_corrections_index", side_effect=RuntimeError("offline")) as load:
            first = corrections.get_corrections_index()
            second = corrections.get_corrections_index()

        self.assertEqual(first, {"by_id": {}, "by_funding": []})
        self.assertIs(first, second)
        load.assert_called_once_with()


if __name__ == "__main__":
    unittest.main()
