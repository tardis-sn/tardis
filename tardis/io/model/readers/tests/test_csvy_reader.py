from collections.abc import Callable
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


@pytest.fixture
def csvy_with_isotopes_fname() -> Path:
    return Path(__file__).parent / "data" / "csvy_with_isotopes.csvy"


@pytest.fixture
def csvy_no_velocity_fname() -> Path:
    return Path(__file__).parent / "data" / "csvy_no_velocity.csvy"


@pytest.fixture
def csvy_no_density_fname() -> Path:
    return Path(__file__).parent / "data" / "csvy_no_density.csvy"


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


@pytest.fixture
def write_csvy(tmp_path: Path) -> Callable:
    def write(text: str) -> Path:
        fname = tmp_path / "model.csvy"
        fname.write_text(text)
        return fname

    return write


def test_load_csvy_converts_velocity_to_cgs_and_drops_inner_boundary_row(
    csvy_with_isotopes_fname,
):
    csvy_data = csvy.load_csvy(csvy_with_isotopes_fname)

    npt.assert_allclose(csvy_data.velocity, [1e9, 1.1e9, 1.2e9])
    assert csvy_data.density is None
    npt.assert_allclose(csvy_data.mass_fractions.loc[1].values, [0.9, 0.7])
    assert list(csvy_data.mass_fractions.columns) == [0, 1]


def test_load_csvy_parses_isotope_mass_fractions(csvy_with_isotopes_fname):
    csvy_data = csvy.load_csvy(csvy_with_isotopes_fname)

    isotope_mass_fractions = csvy_data.isotope_mass_fractions
    assert list(isotope_mass_fractions.index) == [(28, 56)]
    npt.assert_allclose(isotope_mass_fractions.loc[(28, 56)].values, [0.1, 0.3])
    assert list(isotope_mass_fractions.columns) == [0, 1]


def test_load_csvy_without_velocity_raises(csvy_no_velocity_fname):
    with pytest.raises(ValueError, match="Velocity information not found"):
        csvy.load_csvy(csvy_no_velocity_fname)


def test_load_csvy_without_closing_delimiter_raises(write_csvy):
    fname = write_csvy("---\nname: unterminated\n")

    with pytest.raises(ValueError, match="not found"):
        csvy.load_csvy(fname)


def test_load_yaml_from_csvy_returns_yaml_metadata(csvy_with_isotopes_fname):
    yaml_dict = csvy.load_yaml_from_csvy(csvy_with_isotopes_fname)

    assert yaml_dict["name"] == "csvy_with_isotopes"
    assert [field["name"] for field in yaml_dict["datatype"]["fields"]] == [
        "velocity",
        "H",
        "Ni56",
    ]


def test_load_yaml_from_csvy_requires_leading_delimiter(write_csvy):
    with pytest.raises(ValueError, match="First line"):
        csvy.load_yaml_from_csvy(write_csvy("name: no_delimiter\n---\n"))


def test_load_yaml_from_csvy_without_closing_delimiter_raises(write_csvy):
    with pytest.raises(ValueError, match="not found"):
        csvy.load_yaml_from_csvy(write_csvy("---\nname: unterminated\n"))


def test_load_csv_from_csvy_returns_csv_table(csvy_with_isotopes_fname):
    data = csvy.load_csv_from_csvy(csvy_with_isotopes_fname)

    assert list(data.columns) == ["velocity", "H", "Ni56"]
    npt.assert_allclose(data["velocity"], [10000, 11000, 12000])
    npt.assert_allclose(data["Ni56"], [0.0, 0.1, 0.3])


def test_load_csv_from_csvy_without_csv_section_returns_none(write_csvy):
    assert (
        csvy.load_csv_from_csvy(write_csvy("---\nname: empty\n---\n")) is None
    )


def test_load_csv_from_csvy_requires_leading_delimiter(write_csvy):
    with pytest.raises(AssertionError, match="First line"):
        csvy.load_csv_from_csvy(write_csvy("name: no_delimiter\n---\na,b\n"))


def test_load_csv_from_csvy_without_closing_delimiter_raises(write_csvy):
    with pytest.raises(ValueError, match="not found"):
        csvy.load_csv_from_csvy(write_csvy("---\nname: unterminated\n"))


def test_load_csvy_without_density_leaves_density_unset(csvy_no_density_fname):
    assert csvy.load_csvy(csvy_no_density_fname).density is None
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
