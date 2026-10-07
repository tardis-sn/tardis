from typing import Any

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest
from tardisbase.testing.regression_data.regression_data import RegressionData

from tardis.plasma.properties.atomic import (
    IonizationData,
    Levels,
)
from tardis.plasma.properties.general import (
    BetaElectron,
    BetaRadiation,
    GElectron,
    ThermalGElectron,
)
from tardis.plasma.properties.ion_number_density import (
    IonNumberDensity,
    PhiSahaLTE,
    SahaFactor,
    ThermalPhiSahaLTE,
)
from tardis.plasma.properties.level_number_density import LevelNumberDensity
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
    levels, excitation_energy, metastability, g = Levels(None).calculate(
        atom_data, selected_atoms
    )

    t_rad = state["t_rad"].to_numpy()
    t_electrons = t_rad * state["link_t_rad_t_electron"]
    beta_rad = BetaRadiation(None).calculate(t_rad)
    beta_electron = BetaElectron(None).calculate(t_electrons)
    g_electron = GElectron(None).calculate(beta_rad)
    thermal_g_electron = ThermalGElectron(None).calculate(beta_electron)

    atomic_mass = atom_data.atom_data.loc[selected_atoms, "mass"]
    number_density = state["abundance"].mul(state["density"], axis=1).div(
        atomic_mass, axis=0
    )
    ionization_data = IonizationData(None).calculate(atom_data, selected_atoms)

    return {
        "levels": levels,
        "excitation_energy": excitation_energy,
        "metastability": metastability,
        "g": g,
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
    regression_data: RegressionData,
) -> None:
    inputs = lte_equilibrium_inputs
    standard_rad_bf = LevelBoltzmannFactorLTE(None).calculate(
        inputs["excitation_energy"],
        inputs["g"],
        inputs["beta_rad"],
        inputs["levels"],
    )
    expected_rad_bf = regression_data.sync_dataframe(
        standard_rad_bf, key="frame_0"
    )
    pdt.assert_frame_equal(
        standard_rad_bf, expected_rad_bf, rtol=1e-12, atol=0.0
    )

    standard_thermal_bf = ThermalLevelBoltzmannFactorLTE(None).calculate(
        inputs["excitation_energy"],
        inputs["g"],
        inputs["beta_electron"],
        inputs["levels"],
    )
    expected_thermal_bf = regression_data.sync_dataframe(
        standard_thermal_bf, key="frame_1"
    )
    pdt.assert_frame_equal(
        standard_thermal_bf, expected_thermal_bf, rtol=1e-12, atol=0.0
    )

    standard_dilute_bf = LevelBoltzmannFactorDiluteLTE(None).calculate(
        inputs["levels"],
        inputs["g"],
        inputs["excitation_energy"],
        inputs["beta_rad"],
        inputs["w"],
        inputs["metastability"],
    )
    expected_dilute_bf = regression_data.sync_dataframe(
        standard_dilute_bf, key="frame_2"
    )
    pdt.assert_frame_equal(
        standard_dilute_bf, expected_dilute_bf, rtol=1e-12, atol=0.0
    )

    standard_partition = PartitionFunction(None).calculate(standard_rad_bf)
    expected_partition = regression_data.sync_dataframe(
        standard_partition, key="frame_3"
    )
    pdt.assert_frame_equal(
        standard_partition, expected_partition, rtol=1e-12, atol=0.0
    )

    standard_thermal_partition = ThermalLTEPartitionFunction(None).calculate(
        standard_thermal_bf
    )
    expected_thermal_partition = regression_data.sync_dataframe(
        standard_thermal_partition, key="frame_4"
    )
    pdt.assert_frame_equal(
        standard_thermal_partition,
        expected_thermal_partition,
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
    regression_data: RegressionData,
) -> None:
    inputs = lte_equilibrium_inputs
    standard_rad_bf = LevelBoltzmannFactorLTE(None).calculate(
        inputs["excitation_energy"],
        inputs["g"],
        inputs["beta_rad"],
        inputs["levels"],
    )
    standard_partition = PartitionFunction(None).calculate(standard_rad_bf)
    standard_phi = PhiSahaLTE(None).calculate(
        inputs["g_electron"],
        inputs["beta_rad"],
        standard_partition,
        inputs["ionization_data"],
    )
    expected_phi = regression_data.sync_dataframe(standard_phi, key="frame_0")
    pdt.assert_frame_equal(standard_phi, expected_phi, rtol=1e-12, atol=0.0)

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
    expected_thermal_phi = regression_data.sync_dataframe(
        standard_thermal_phi, key="frame_1"
    )
    pdt.assert_frame_equal(
        standard_thermal_phi, expected_thermal_phi, rtol=1e-12, atol=0.0
    )

    standard_phi_ik = SahaFactor(None).calculate(
        standard_thermal_phi, standard_thermal_bf, standard_thermal_partition
    )
    expected_phi_ik = regression_data.sync_dataframe(
        standard_phi_ik, key="frame_2"
    )
    pdt.assert_frame_equal(
        standard_phi_ik, expected_phi_ik, rtol=1e-12, atol=0.0
    )


def test_lte_ion_and_level_populations_conserve_elements(
    lte_equilibrium_inputs: dict[str, Any],
    regression_data: RegressionData,
) -> None:
    inputs = lte_equilibrium_inputs
    level_bf = LevelBoltzmannFactorLTE(None).calculate(
        inputs["excitation_energy"],
        inputs["g"],
        inputs["beta_rad"],
        inputs["levels"],
    )
    partition = PartitionFunction(None).calculate(level_bf)
    phi = PhiSahaLTE(None).calculate(
        inputs["g_electron"],
        inputs["beta_rad"],
        partition,
        inputs["ionization_data"],
    )

    standard_ions, standard_electrons = IonNumberDensity(None).calculate(
        phi, partition, inputs["number_density"]
    )
    expected_ions = regression_data.sync_dataframe(standard_ions, key="frame_0")
    pdt.assert_frame_equal(
        standard_ions, expected_ions, rtol=1e-12, atol=0.0
    )
    standard_electron_frame = pd.DataFrame({"value": standard_electrons})
    expected_electrons = regression_data.sync_dataframe(
        standard_electron_frame, key="allclose_0"
    )
    npt.assert_allclose(
        standard_electron_frame.to_numpy(),
        expected_electrons.to_numpy(),
        rtol=5e-2,
    )

    standard_levels = LevelNumberDensity(None).calculate(
        level_bf, standard_ions, inputs["levels"], partition
    )
    expected_levels = regression_data.sync_dataframe(standard_levels, key="frame_1")
    pdt.assert_frame_equal(
        standard_levels, expected_levels, rtol=1e-12, atol=0.0
    )

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
    # Level number densities are normalized by their partition function; this
    # identity is also independent of the ion-density solver tolerance.
    npt.assert_allclose(
        level_by_ion.to_numpy(), standard_ions.to_numpy(), rtol=1e-12
    )
    assert (standard_ions >= 0).all().all()
    assert (standard_levels >= 0).all().all()
