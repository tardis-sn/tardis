from pathlib import Path
import tardis
from tardis.io.model.readers import csvy
from tardis.io.configuration.config_validator import validate_dict
from jsonschema import exceptions as json_schema_exc
import pytest
import os
from astropy import units as u
import numpy.testing as npt


@pytest.fixture
def csvy_full_fname(example_model_file_dir: Path):
    return example_model_file_dir / "csvy_full.csvy"


@pytest.fixture
def csvy_nocsv_fname(example_model_file_dir: Path):
    return example_model_file_dir / "csvy_nocsv.csvy"


@pytest.fixture
def csvy_missing_fname(example_model_file_dir: Path):
    return example_model_file_dir / "csvy_missing.csvy"


def test_csvy_finds_csv_first_line(csvy_full_fname):
    yaml_dict, csv = csvy.load_csvy(csvy_full_fname)
    npt.assert_almost_equal(csv["velocity"][0], 10000)


def test_csv_colnames_equiv_datatype_fields(csvy_full_fname):
    yaml_dict, csv = csvy.load_csvy(csvy_full_fname)
    datatype_names = [od["name"] for od in yaml_dict["datatype"]["fields"]]
    for key in csv.columns:
        assert key in datatype_names
    for name in datatype_names:
        assert name in csv.columns


def test_csvy_nocsv_data_is_none(csvy_nocsv_fname):
    yaml_dict, csv = csvy.load_csvy(csvy_nocsv_fname)
    assert csv is None


def test_missing_required_property(csvy_missing_fname):
    yaml_dict, csv = csvy.load_csvy(csvy_missing_fname)
    with pytest.raises(Exception):
        vy = validate_dict(
            yaml_dict,
            schemapath=os.path.join(
                tardis.__path__[0], "io", "schemas", "csvy_model.yml"
            ),
        )
