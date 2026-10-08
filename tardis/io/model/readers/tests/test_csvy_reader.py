from pathlib import Path

import numpy.testing as npt
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


CSVY_WITH_ISOTOPES = """---
name: csvy_with_isotopes
model_density_time_0: 1 day
model_isotope_time_0: 0 day
description: Velocity and abundances (including an isotope) in the CSV section.
tardis_model_config_version: v1.0
datatype:
  fields:
    - name: velocity
      unit: km/s
      desc: velocities of shell outer boundaries.
    - name: H
      desc: fractional H abundance
    - name: Ni56
      desc: fractional Ni56 abundance

density:
  type: uniform
  value: 1e-10 g/cm^3
---
velocity,H,Ni56
10000,0.0,0.0
11000,0.9,0.1
12000,0.7,0.3
"""


@pytest.fixture
def write_csvy(tmp_path: Path):
    def write(text: str) -> Path:
        fname = tmp_path / "model.csvy"
        fname.write_text(text)
        return fname

    return write


def test_load_csvy_converts_velocity_to_cgs_and_drops_inner_boundary_row(
    write_csvy,
):
    csvy_data = csvy.load_csvy(write_csvy(CSVY_WITH_ISOTOPES))

    npt.assert_allclose(csvy_data.velocity, [1e9, 1.1e9, 1.2e9])
    assert csvy_data.density is None
    npt.assert_allclose(csvy_data.mass_fractions.loc[1].values, [0.9, 0.7])
    assert list(csvy_data.mass_fractions.columns) == [0, 1]


def test_load_csvy_parses_isotope_mass_fractions(write_csvy):
    csvy_data = csvy.load_csvy(write_csvy(CSVY_WITH_ISOTOPES))

    isotope_mass_fractions = csvy_data.isotope_mass_fractions
    assert list(isotope_mass_fractions.index) == [(28, 56)]
    npt.assert_allclose(isotope_mass_fractions.loc[(28, 56)].values, [0.1, 0.3])
    assert list(isotope_mass_fractions.columns) == [0, 1]


def test_load_csvy_without_velocity_raises(write_csvy):
    fname = write_csvy(
        """---
name: no_velocity
model_density_time_0: 1 day
model_isotope_time_0: 0 day
description: No velocity information anywhere.
tardis_model_config_version: v1.0
density:
  type: uniform
  value: 1e-10 g/cm^3
---
"""
    )

    with pytest.raises(ValueError, match="Velocity information not found"):
        csvy.load_csvy(fname)


def test_load_csvy_without_closing_delimiter_raises(write_csvy):
    fname = write_csvy("---\nname: unterminated\n")

    with pytest.raises(ValueError, match="not found"):
        csvy.load_csvy(fname)


def test_load_yaml_from_csvy_returns_yaml_metadata(write_csvy):
    yaml_dict = csvy.load_yaml_from_csvy(write_csvy(CSVY_WITH_ISOTOPES))

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


def test_load_csv_from_csvy_returns_csv_table(write_csvy):
    data = csvy.load_csv_from_csvy(write_csvy(CSVY_WITH_ISOTOPES))

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


def test_load_csvy_without_density_leaves_density_unset(write_csvy):
    fname = write_csvy(
        """---
name: no_density
model_density_time_0: 1 day
model_isotope_time_0: 0 day
description: Velocity and abundances only.
tardis_model_config_version: v1.0
datatype:
  fields:
    - name: velocity
      unit: km/s
      desc: velocities of shell boundaries.
    - name: H
      desc: fractional H abundance
---
velocity,H
10000,0.0
11000,1.0
"""
    )

    assert csvy.load_csvy(fname).density is None
