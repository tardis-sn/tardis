"""
Sanity check the tests in the suite.

To run this, unless you otherwise know what you're doing:

    * install ``pipx`` (see its documentation page)
    * install nox via ``pipx install nox``
    * run ``nox`` in the root of this repository
"""

from pathlib import Path
import json

from jsonschema.validators import validator_for
import pytest

ROOT = Path(__file__).parent
VERSIONS = ROOT / "tests"
SPECS = json.loads(VERSIONS.joinpath("specifications.json").read_text())

_SCHEMA = json.loads(ROOT.joinpath("test-schema.json").read_text())
Validator = validator_for(_SCHEMA)
Validator.check_schema(_SCHEMA)
VALIDATOR = Validator(_SCHEMA)


@pytest.mark.parametrize("test_path", VERSIONS.glob("*/**/*.json"))
def test_tests_are_valid(test_path):
    try:
        test = json.loads(test_path.read_text())
    except json.JSONDecodeError:
        pytest.fail(f"{test_path} contains invalid JSON")
    else:
        VALIDATOR.validate(test)


@pytest.mark.parametrize(
    "version",
    [version for version in VERSIONS.iterdir() if version.is_dir()],
)
def test_specification_directories_are_identified(version):
    assert version.name in SPECS
