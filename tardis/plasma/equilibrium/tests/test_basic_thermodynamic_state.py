from typing import Any

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
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
from tardis.iip_plasma.properties.plasma_input import JBlues as IIPJBlues
from tardis.plasma.properties.atomic import IonizationData, Levels
from tardis.plasma.properties.general import (
    BetaElectron,
    BetaRadiation,
    ElectronTemperature,
    GElectron,
)
from tardis.plasma.properties.plasma_input import JBlues as StandardJBlues
from tardis.util.base import intensity_black_body



def test_atomic_data_basic_properties_match_iip(
    basic_thermodynamic_state: dict[str, Any],
) -> None:
    state = basic_thermodynamic_state
    atom_data = state["atomic_data"]
    selected_atoms = state["selected_atoms"]

    standard_levels = Levels(None).calculate(atom_data, selected_atoms)
    iip_levels = IIPLevels(None).calculate(atom_data, selected_atoms)
    pdt.assert_index_equal(standard_levels[0], iip_levels[0])
    for standard, iip in zip(standard_levels[1:], iip_levels[1:], strict=True):
        pdt.assert_series_equal(standard, iip)

    standard_ionization = IonizationData(None).calculate(
        atom_data, selected_atoms
    )
    iip_ionization = IIPIonizationData(None).calculate(
        atom_data, selected_atoms
    )
    pd.testing.assert_series_equal(standard_ionization, iip_ionization)

    standard_mass = atom_data.atom_data.loc[selected_atoms, "mass"]
    iip_mass = IIPAtomicMass(None).calculate(atom_data, selected_atoms)
    pd.testing.assert_series_equal(standard_mass, iip_mass)


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


def test_thermodynamic_inputs_match_iip_and_satisfy_identities(
    basic_thermodynamic_state: dict[str, Any],
) -> None:
    state = basic_thermodynamic_state
    t_rad = state["t_rad"].to_numpy()
    link = state["link_t_rad_t_electron"]

    standard_t_electrons = ElectronTemperature(None).calculate(t_rad, link)
    iip_t_electrons = IIPElectronTemperature(None).calculate(t_rad, link)
    npt.assert_allclose(standard_t_electrons, iip_t_electrons)
    npt.assert_allclose(standard_t_electrons, link * t_rad)

    standard_beta_rad = BetaRadiation(None).calculate(t_rad)
    iip_beta_rad = IIPBetaRadiation(None).calculate(t_rad)
    standard_beta_electron = BetaElectron(None).calculate(standard_t_electrons)
    iip_beta_electron = IIPBetaElectron(None).calculate(iip_t_electrons)
    npt.assert_allclose(standard_beta_rad, iip_beta_rad, rtol=3e-7)
    npt.assert_allclose(standard_beta_electron, iip_beta_electron, rtol=3e-7)

    standard_g = GElectron(None).calculate(standard_beta_rad)
    iip_g = IIPGElectron(None).calculate(iip_beta_rad)
    expected_g = (
        2
        * np.pi
        * const.m_e.cgs.value
        / standard_beta_rad
        / const.h.cgs.value**2
    ) ** 1.5
    npt.assert_allclose(standard_g, iip_g, rtol=5e-7)
    npt.assert_allclose(standard_g, expected_g)


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


def test_line_mean_intensity_inputs_preserve_canonical_shell_order(
    basic_thermodynamic_state: dict[str, Any],
) -> None:
    state = basic_thermodynamic_state
    lines = (
        state["atomic_data"]
        .lines.loc[(1, 0, slice(None), slice(None))]
        .iloc[:2]
    )
    line_index = lines.index
    expected = state["radiation_field"].calculate_mean_intensity(
        lines.nu.to_numpy() * u.Hz
    )
    supplied = pd.DataFrame(
        expected,
        index=line_index,
        columns=pd.Index([0, 1, 2], name="shell"),
    )

    standard_input = StandardJBlues()
    standard_input.set_value(supplied)
    iip_input = IIPJBlues()
    iip_input.set_value(supplied)

    pdt.assert_frame_equal(standard_input.j_blues, supplied)
    npt.assert_array_equal(iip_input.j_blues, standard_input.j_blues.to_numpy())
