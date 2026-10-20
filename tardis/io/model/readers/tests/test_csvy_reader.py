from pathlib import Path

import numpy.testing as npt
import pandas as pd
import pytest

from tardis.io.model import csvy


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
    csvy_data = csvy.load_csvy(csvy_full_fname)

    npt.assert_almost_equal(csvy_data.raw_csv_data["velocity"][0], 10000)


def test_csv_colnames_equiv_datatype_fields(csvy_full_fname):
    csvy_data = csvy.load_csvy(csvy_full_fname)

    datatype_names = [
        od["name"] for od in csvy_data.model_config.datatype.fields
    ]
    for key in csvy_data.raw_csv_data.columns:
        assert key in datatype_names
    for name in datatype_names:
        assert name in csvy_data.raw_csv_data.columns


def test_csvy_nocsv_data_is_none(csvy_nocsv_fname):
    csvy_data = csvy.load_csvy(csvy_nocsv_fname)
    assert csvy_data.raw_csv_data is None


def test_missing_required_property(csvy_missing_fname):
    # Validation now happens inside load_csvy, so it should raise during loading
    with pytest.raises(Exception):
        csvy_data = csvy.load_csvy(csvy_missing_fname)


def test_parse_csv_mass_fractions_ignores_capitalization():
    csv_data = pd.DataFrame(
        {
            "velocity": [9000.0, 10000.0],
            "h": [0.5, 0.3],
            "HE": [0.3, 0.3],
            "Si": [0.1, 0.2],
            "ni56": [0.1, 0.2],
        }
    )
    _, mass_fractions, isotope_mass_fractions = csvy.parse_csv_mass_fractions(
        csv_data
    )
    npt.assert_allclose(mass_fractions.loc[1].values, [0.5, 0.3])
    npt.assert_allclose(mass_fractions.loc[2].values, [0.3, 0.3])
    npt.assert_allclose(mass_fractions.loc[14].values, [0.1, 0.2])
    npt.assert_allclose(isotope_mass_fractions.loc[(28, 56)].values, [0.1, 0.2])
