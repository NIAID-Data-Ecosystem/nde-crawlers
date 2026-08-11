import importlib.util
import logging
import sys
import types
import unittest
from pathlib import Path


HUB_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = HUB_ROOT / "utils" / "validate.py"
sys.path.insert(0, str(HUB_ROOT))

config = sys.modules.get("config")
if config is None:
    config = types.ModuleType("config")
    sys.modules["config"] = config
config.logger = logging.getLogger("test_validate")

package_name = "validate_test_utils"
package = types.ModuleType(package_name)
package.__path__ = [str(HUB_ROOT / "utils")]
sys.modules[package_name] = package

common = types.ModuleType(f"{package_name}.common")
common.as_list = lambda value: [] if value is None else value if isinstance(value, list) else [value]
sys.modules[common.__name__] = common

spec = importlib.util.spec_from_file_location(f"{package_name}.validate", MODULE_PATH)
validate = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = validate
spec.loader.exec_module(validate)


class SchemaValidationTests(unittest.TestCase):
    def test_archived_at_failure_names_the_record(self):
        doc = {
            "_id": "PMC123",
            "@type": "Dataset",
            "url": "https://example.org/PMC123",
            "includedInDataCatalog": {"name": "Example", "archivedAt": None},
        }

        with self.assertRaisesRegex(AssertionError, r"archivedAt.*record _id='PMC123'"):
            validate.check_schema(doc)


if __name__ == "__main__":
    unittest.main()
