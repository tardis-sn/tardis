from typing import Any

import numpy as np
import numpy.testing as npt
import pandas.testing as pdt
import pytest

from tardis.iip_plasma.properties.atomic import AtomicMass, Levels
from tardis.iip_plasma.properties.general import (
    NumberDensity,
)
from tardis.iip_plasma.properties.ion_population import (
    IonNumberDensity as IIPIonNumberDensity,
)
from tardis.iip_plasma.properties.ion_population import (
    PhiSahaElectrons,
    PhiSahaLTE,
)
from tardis.iip_plasma.properties.level_population import (
    LevelNumberDensity as IIPLevelNumberDensity,
)
from tardis.iip_plasma.properties.level_population import PhiLucy
from tardis.iip_plasma.properties.partition_function import (
    LevelBoltzmannFactorDiluteLTE as IIPLevelBoltzmannFactorDiluteLTE,
)
from tardis.iip_plasma.properties.partition_function import (
    LevelBoltzmannFactorLTE as IIPLevelBoltzmannFactorLTE,
)
from tardis.iip_plasma.properties.partition_function import (
    LevelBoltzmannFactorLTETe,
)
from tardis.iip_plasma.properties.partition_function import (
    PartitionFunction as IIPPartitionFunction,
)
from tardis.plasma.properties.atomic import (
    IonizationData as StandardIonizationData,
)
from tardis.plasma.properties.atomic import Levels as StandardLevels
from tardis.plasma.properties.general import (
    BetaElectron as StandardBetaElectron,
)
from tardis.plasma.properties.general import (
    BetaRadiation as StandardBetaRadiation,
)
from tardis.plasma.properties.general import (
    GElectron as StandardGElectron,
)
from tardis.plasma.properties.general import (
    ThermalGElectron,
)
from tardis.plasma.properties.ion_population import (
    IonNumberDensity as StandardIonNumberDensity,
)
from tardis.plasma.properties.ion_population import (
    PhiSahaLTE as StandardPhiSahaLTE,
)
from tardis.plasma.properties.ion_population import (
    SahaFactor,
    ThermalPhiSahaLTE,
)
from tardis.plasma.properties.level_population import (
    LevelNumberDensity as StandardLevelNumberDensity,
)
from tardis.plasma.properties.partition_function import (
    LevelBoltzmannFactorDiluteLTE,
    LevelBoltzmannFactorLTE,
    PartitionFunction,
    ThermalLevelBoltzmannFactorLTE,
    ThermalLTEPartitionFunction,
)


@pytest.fixture
def lte_equilibrium_inputs(
    basic_thermodynamic_state: dict[str, Any],
) -> dict[str, Any]:
    state = basic_thermodynamic_state
    atom_data = state["atomic_data"]
    selected_atoms = state["selected_atoms"]
    levels, excitation_energy, metastability, g = StandardLevels(
        None
    ).calculate(atom_data, selected_atoms)
    iip_levels, iip_excitation_energy, iip_metastability, iip_g = Levels(
        None
    ).calculate(atom_data, selected_atoms)

    t_rad = state["t_rad"].to_numpy()
    t_electrons = t_rad * state["link_t_rad_t_electron"]
    beta_rad = StandardBetaRadiation(None).calculate(t_rad)
    beta_electron = StandardBetaElectron(None).calculate(t_electrons)
    g_electron = StandardGElectron(None).calculate(beta_rad)
    thermal_g_electron = ThermalGElectron(None).calculate(beta_electron)

    masses = AtomicMass(None).calculate(atom_data, selected_atoms)
    number_density = NumberDensity(None).calculate(
        masses, state["abundance"], state["density"]
    )
    ionization_data = StandardIonizationData(None).calculate(
        atom_data, selected_atoms
    )

    return {
        "levels": levels,
        "excitation_energy": excitation_energy,
        "metastability": metastability,
        "g": g,
        "iip_levels": iip_levels,
        "iip_excitation_energy": iip_excitation_energy,
        "iip_metastability": iip_metastability,
        "iip_g": iip_g,
        "beta_rad": beta_rad,
        "beta_electron": beta_electron,
        "g_electron": g_electron,
        "thermal_g_electron": thermal_g_electron,
        "number_density": number_density,
        "ionization_data": ionization_data,
        "w": state["dilution_factor"].to_numpy(),
    }


def test_boltzmann_factors_and_partition_functions_match_iip(
    lte_equilibrium_inputs: dict[str, Any],
) -> None:
    inputs = lte_equilibrium_inputs

    standard_rad_bf = LevelBoltzmannFactorLTE(None).calculate(
        inputs["excitation_energy"],
        inputs["g"],
        inputs["beta_rad"],
        inputs["levels"],
    )
    iip_rad_bf = IIPLevelBoltzmannFactorLTE(None).calculate(
        inputs["iip_excitation_energy"],
        inputs["iip_g"],
        inputs["beta_rad"],
        inputs["iip_levels"],
    )
    pdt.assert_frame_equal(standard_rad_bf, iip_rad_bf, rtol=1e-12, atol=0.0)

    standard_thermal_bf = ThermalLevelBoltzmannFactorLTE(None).calculate(
        inputs["excitation_energy"],
        inputs["g"],
        inputs["beta_electron"],
        inputs["levels"],
    )
    iip_thermal_bf = LevelBoltzmannFactorLTETe(None).calculate(
        inputs["iip_excitation_energy"],
        inputs["iip_g"],
        inputs["beta_electron"],
        inputs["iip_levels"],
    )
    pdt.assert_frame_equal(
        standard_thermal_bf, iip_thermal_bf, rtol=1e-12, atol=0.0
    )

    standard_dilute_bf = LevelBoltzmannFactorDiluteLTE(None).calculate(
        inputs["levels"],
        inputs["g"],
        inputs["excitation_energy"],
        inputs["beta_rad"],
        inputs["w"],
        inputs["metastability"],
    )
    iip_dilute_bf = IIPLevelBoltzmannFactorDiluteLTE(None).calculate(
        inputs["iip_levels"],
        inputs["iip_g"],
        inputs["iip_excitation_energy"],
        inputs["beta_rad"],
        inputs["w"],
        inputs["iip_metastability"],
    )
    pdt.assert_frame_equal(
        standard_dilute_bf, iip_dilute_bf, rtol=1e-12, atol=0.0
    )

    standard_partition = PartitionFunction(None).calculate(standard_rad_bf)
    iip_partition = IIPPartitionFunction(None).calculate(iip_rad_bf)
    pdt.assert_frame_equal(
        standard_partition, iip_partition, rtol=1e-12, atol=0.0
    )

    standard_thermal_partition = ThermalLTEPartitionFunction(None).calculate(
        standard_thermal_bf
    )
    iip_thermal_partition = IIPPartitionFunction(None).calculate(iip_thermal_bf)
    pdt.assert_frame_equal(
        standard_thermal_partition,
        iip_thermal_partition,
        rtol=1e-12,
        atol=0.0,
    )

    expected_partition = standard_rad_bf.groupby(
        level=["atomic_number", "ion_number"]
    ).sum()
    pdt.assert_frame_equal(standard_partition, expected_partition)
    assert (standard_partition > 0).all().all()


def test_dilute_lte_correction_only_changes_non_metastable_levels(
    lte_equilibrium_inputs: dict[str, Any],
) -> None:
    inputs = lte_equilibrium_inputs
    dilute = LevelBoltzmannFactorDiluteLTE(None).calculate(
        inputs["levels"],
        inputs["g"],
        inputs["excitation_energy"],
        inputs["beta_rad"],
        inputs["w"],
        inputs["metastability"],
    )
    ordinary = LevelBoltzmannFactorLTE(None).calculate(
        inputs["excitation_energy"],
        inputs["g"],
        inputs["beta_rad"],
        inputs["levels"],
    )
    metastable = inputs["metastability"].to_numpy()
    ratio = dilute.to_numpy() / ordinary.to_numpy()
    npt.assert_allclose(ratio[metastable], 1.0)
    npt.assert_allclose(
        ratio[~metastable],
        np.broadcast_to(inputs["w"], ratio[~metastable].shape),
    )


def test_saha_factors_and_phi_ik_match_iip(
    lte_equilibrium_inputs: dict[str, Any],
) -> None:
    inputs = lte_equilibrium_inputs
    standard_rad_bf = LevelBoltzmannFactorLTE(None).calculate(
        inputs["excitation_energy"],
        inputs["g"],
        inputs["beta_rad"],
        inputs["levels"],
    )
    standard_partition = PartitionFunction(None).calculate(standard_rad_bf)
    standard_phi = StandardPhiSahaLTE(None).calculate(
        inputs["g_electron"],
        inputs["beta_rad"],
        standard_partition,
        inputs["ionization_data"],
    )
    iip_phi = PhiSahaLTE(None).calculate(
        inputs["g_electron"],
        inputs["beta_rad"],
        standard_partition,
        inputs["ionization_data"],
    )
    pdt.assert_frame_equal(standard_phi, iip_phi, rtol=1e-12, atol=0.0)

    standard_thermal_bf = ThermalLevelBoltzmannFactorLTE(None).calculate(
        inputs["excitation_energy"],
        inputs["g"],
        inputs["beta_electron"],
        inputs["levels"],
    )
    standard_thermal_partition = ThermalLTEPartitionFunction(None).calculate(
        standard_thermal_bf
    )
    standard_thermal_phi = ThermalPhiSahaLTE(None).calculate(
        inputs["thermal_g_electron"],
        inputs["beta_electron"],
        standard_thermal_partition,
        inputs["ionization_data"],
    )
    iip_thermal_phi = PhiSahaElectrons(None).calculate(
        inputs["thermal_g_electron"],
        inputs["beta_electron"],
        standard_thermal_partition,
        inputs["ionization_data"],
    )
    pdt.assert_frame_equal(
        standard_thermal_phi, iip_thermal_phi, rtol=1e-12, atol=0.0
    )

    standard_phi_ik = SahaFactor(None).calculate(
        standard_thermal_phi, standard_thermal_bf, standard_thermal_partition
    )
    iip_phi_lucy = PhiLucy(None).calculate(
        iip_thermal_phi, standard_thermal_bf, standard_thermal_partition
    )
    pdt.assert_frame_equal(standard_phi_ik, iip_phi_lucy, rtol=1e-12, atol=0.0)


def test_lte_ion_and_level_populations_conserve_elements(
    lte_equilibrium_inputs: dict[str, Any],
) -> None:
    inputs = lte_equilibrium_inputs
    level_bf = LevelBoltzmannFactorLTE(None).calculate(
        inputs["excitation_energy"],
        inputs["g"],
        inputs["beta_rad"],
        inputs["levels"],
    )
    partition = PartitionFunction(None).calculate(level_bf)
    phi = StandardPhiSahaLTE(None).calculate(
        inputs["g_electron"],
        inputs["beta_rad"],
        partition,
        inputs["ionization_data"],
    )

    standard_ions, standard_electrons = StandardIonNumberDensity(
        None
    ).calculate(phi, partition, inputs["number_density"])
    iip_ions, iip_electrons = IIPIonNumberDensity(None).calculate(
        phi, partition, inputs["number_density"]
    )
    pdt.assert_frame_equal(standard_ions, iip_ions, rtol=1e-12, atol=0.0)
    # Both legacy ion solvers stop when their electron-density iteration
    # changes by less than 5%; this is a solver-parity check, not a
    # conservation tolerance.
    npt.assert_allclose(standard_electrons, iip_electrons, rtol=5e-2)

    standard_levels = StandardLevelNumberDensity(None).calculate(
        level_bf, standard_ions, inputs["levels"], partition
    )
    iip_levels = IIPLevelNumberDensity(None).calculate(
        level_bf, iip_ions, inputs["levels"], partition
    )
    pdt.assert_frame_equal(standard_levels, iip_levels, rtol=1e-12, atol=0.0)

    ion_by_element = standard_ions.groupby(level="atomic_number").sum()
    pdt.assert_index_equal(ion_by_element.index, inputs["number_density"].index)
    pdt.assert_index_equal(
        ion_by_element.columns,
        inputs["number_density"].columns,
        check_names=False,
    )
    # Element conservation is an algebraic invariant, so only summation
    # round-off is allowed here, independently of solver convergence.
    npt.assert_allclose(
        ion_by_element.to_numpy(),
        inputs["number_density"].to_numpy(),
        rtol=1e-12,
    )
    level_by_ion = standard_levels.groupby(
        level=["atomic_number", "ion_number"]
    ).sum()
    # Level populations are normalized by their partition function; this
    # identity is also independent of the ion-density solver tolerance.
    npt.assert_allclose(
        level_by_ion.to_numpy(), standard_ions.to_numpy(), rtol=1e-12
    )
    assert (standard_ions >= 0).all().all()
    assert (standard_levels >= 0).all().all()
