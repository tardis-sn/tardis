from typing import Any

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest
from astropy import units as u

from tardis import constants as const
from tardis.iip_plasma.properties.atomic import (
    AtomicMass as IIPAtomicMass,
)
from tardis.iip_plasma.properties.atomic import (
    IonizationData as IIPIonizationData,
)
from tardis.iip_plasma.properties.atomic import (
    Levels as IIPLevels,
)
from tardis.iip_plasma.properties.general import (
    BetaElectron as IIPBetaElectron,
)
from tardis.iip_plasma.properties.general import (
    BetaRadiation as IIPBetaRadiation,
)
from tardis.iip_plasma.properties.general import (
    ElectronTemperature as IIPElectronTemperature,
)
from tardis.iip_plasma.properties.general import (
    GElectron as IIPGElectron,
)
from tardis.iip_plasma.properties.general import (
    NumberDensity as IIPNumberDensity,
)
from tardis.plasma.properties.atomic import IonizationData, Levels
from tardis.plasma.properties.general import (
    BetaElectron,
    BetaRadiation,
    ElectronTemperature,
    GElectron,
)
from tardis.util.base import intensity_black_body


@pytest.fixture
def atomic_property_values(
    basic_thermodynamic_state: dict[str, Any],
) -> dict[str, Any]:
    state = basic_thermodynamic_state
    atom_data = state["atomic_data"]
    selected_atoms = state["selected_atoms"]

    standard_levels = Levels(None).calculate(atom_data, selected_atoms)
    iip_levels = IIPLevels(None).calculate(atom_data, selected_atoms)
    standard_ionization = IonizationData(None).calculate(
        atom_data, selected_atoms
    )
    iip_ionization = IIPIonizationData(None).calculate(
        atom_data, selected_atoms
    )
    standard_mass = atom_data.atom_data.loc[selected_atoms, "mass"]
    iip_mass = IIPAtomicMass(None).calculate(atom_data, selected_atoms)
    return {
        "standard_levels": standard_levels,
        "iip_levels": iip_levels,
        "standard_ionization": standard_ionization,
        "iip_ionization": iip_ionization,
        "standard_mass": standard_mass,
        "iip_mass": iip_mass,
    }


def test_atomic_levels_match_iip(atomic_property_values: dict[str, Any]) -> None:
    pdt.assert_index_equal(
        atomic_property_values["standard_levels"][0],
        atomic_property_values["iip_levels"][0],
    )
    for standard, iip in zip(
        atomic_property_values["standard_levels"][1:],
        atomic_property_values["iip_levels"][1:],
        strict=True,
    ):
        pdt.assert_series_equal(standard, iip)


def test_ionization_data_matches_iip(
    atomic_property_values: dict[str, Any],
) -> None:
    pd.testing.assert_series_equal(
        atomic_property_values["standard_ionization"],
        atomic_property_values["iip_ionization"],
    )


def test_atomic_mass_matches_iip(atomic_property_values: dict[str, Any]) -> None:
    pd.testing.assert_series_equal(
        atomic_property_values["standard_mass"],
        atomic_property_values["iip_mass"],
    )


def test_number_density_and_mass_reconstruct_density(
    basic_thermodynamic_state: dict[str, Any],
) -> None:
    state = basic_thermodynamic_state
    atom_data = state["atomic_data"]
    abundance = state["abundance"]
    density = state["density"]
    masses = IIPAtomicMass(None).calculate(atom_data, state["selected_atoms"])

    number_density = IIPNumberDensity(None).calculate(
        masses, abundance, density
    )
    reconstructed_density = number_density.mul(masses, axis=0).sum(axis=0)

    npt.assert_allclose(reconstructed_density.to_numpy(), density.to_numpy())
    assert (number_density >= 0).all().all()
    assert number_density.index.equals(abundance.index)
    assert number_density.columns.equals(abundance.columns)


@pytest.fixture
def thermodynamic_property_values(
    basic_thermodynamic_state: dict[str, Any],
) -> dict[str, Any]:
    state = basic_thermodynamic_state
    t_rad = state["t_rad"].to_numpy()
    link = state["link_t_rad_t_electron"]

    standard_t_electrons = ElectronTemperature(None).calculate(t_rad, link)
    iip_t_electrons = IIPElectronTemperature(None).calculate(t_rad, link)
    standard_beta_rad = BetaRadiation(None).calculate(t_rad)
    iip_beta_rad = IIPBetaRadiation(None).calculate(t_rad)
    standard_beta_electron = BetaElectron(None).calculate(standard_t_electrons)
    iip_beta_electron = IIPBetaElectron(None).calculate(iip_t_electrons)
    standard_g = GElectron(None).calculate(standard_beta_rad)
    iip_g = IIPGElectron(None).calculate(iip_beta_rad)
    return {
        "t_rad": t_rad,
        "link": link,
        "standard_t_electrons": standard_t_electrons,
        "iip_t_electrons": iip_t_electrons,
        "standard_beta_rad": standard_beta_rad,
        "iip_beta_rad": iip_beta_rad,
        "standard_beta_electron": standard_beta_electron,
        "iip_beta_electron": iip_beta_electron,
        "standard_g": standard_g,
        "iip_g": iip_g,
    }


def test_electron_temperature_matches_iip(
    thermodynamic_property_values: dict[str, Any],
) -> None:
    npt.assert_allclose(
        thermodynamic_property_values["standard_t_electrons"],
        thermodynamic_property_values["iip_t_electrons"],
    )
    npt.assert_allclose(
        thermodynamic_property_values["standard_t_electrons"],
        thermodynamic_property_values["link"]
        * thermodynamic_property_values["t_rad"],
    )


def test_radiation_beta_matches_iip(
    thermodynamic_property_values: dict[str, Any],
) -> None:
    npt.assert_allclose(
        thermodynamic_property_values["standard_beta_rad"],
        thermodynamic_property_values["iip_beta_rad"],
        rtol=3e-7,
    )


def test_electron_beta_matches_iip(
    thermodynamic_property_values: dict[str, Any],
) -> None:
    npt.assert_allclose(
        thermodynamic_property_values["standard_beta_electron"],
        thermodynamic_property_values["iip_beta_electron"],
        rtol=3e-7,
    )


def test_electron_statistical_factor_matches_iip(
    thermodynamic_property_values: dict[str, Any],
) -> None:
    standard_beta_rad = thermodynamic_property_values["standard_beta_rad"]
    expected_g = (
        2 * np.pi * const.m_e.cgs.value / standard_beta_rad / const.h.cgs.value**2
    ) ** 1.5
    npt.assert_allclose(
        thermodynamic_property_values["standard_g"],
        thermodynamic_property_values["iip_g"],
        rtol=5e-7, #  iip_plasma uses raw astropy constants not tardis.constants
    )
    npt.assert_allclose(thermodynamic_property_values["standard_g"], expected_g)


def test_dilute_planckian_mean_intensity_matches_analytic_planck_function(
    basic_thermodynamic_state: dict[str, Any],
) -> None:
    state = basic_thermodynamic_state
    frequencies = (
        state["atomic_data"]
        .lines.loc[(1, 0, slice(None), slice(None)), "nu"]
        .iloc[:3]
        .to_numpy()
        * u.Hz
    )
    actual = state["radiation_field"].calculate_mean_intensity(frequencies)
    expected = state["dilution_factor"].to_numpy() * intensity_black_body(
        frequencies[np.newaxis].T, state["t_rad"].to_numpy() * u.K
    )

    # ``intensity_black_body`` and the radiation-field API return cgs values
    # without an Astropy unit wrapper.
    assert isinstance(actual, np.ndarray)
    npt.assert_allclose(actual, expected)
