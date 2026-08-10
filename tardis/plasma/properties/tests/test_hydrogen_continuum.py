from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.plasma.assembly.base import IIPPlasmaSolverFactory
from tardis.plasma.properties.hydrogen_continuum import (
    HydrogenContinuumLTEPopulations,
    IIPContinuumPopulations,
)
from tardis.plasma.properties.iip_property_collections import (
    hydrogen_continuum_properties,
)
from tardis.plasma.properties.radiative_properties import (
    StimulatedEmissionFactor,
)


def test_iip_continuum_outputs_have_one_coupled_owner() -> None:
    owners = [
        property_class
        for property_class in hydrogen_continuum_properties
        if {
            "hydrogen_continuum_level_boltzmann_factor",
            "hydrogen_continuum_partition_function",
            "ion_number_density",
            "electron_densities",
            "level_number_density",
        }
        & set(property_class.outputs)
    ]

    assert owners == [IIPContinuumPopulations]
    assert "_build_trial_state" not in IIPContinuumPopulations.__dict__


def test_lte_populations_use_current_electron_densities() -> None:
    columns = pd.Index([0], name="shell")
    ion_index = pd.MultiIndex.from_tuples(
        [(1, 0), (1, 1)], names=["atomic_number", "ion_number"]
    )
    level_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (1, 1, 0)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    current_electron_densities = pd.Series([5.0], index=columns)
    number_density = pd.DataFrame(
        [[12.0]], index=pd.Index([1], name="atomic_number"), columns=columns
    )
    partition_function = pd.DataFrame(
        [[1.0], [1.0]], index=ion_index, columns=columns
    )
    thermal_level_factor = pd.DataFrame(
        [[2.0], [1.0]], index=level_index, columns=columns
    )
    ionization_data = pd.Series([0.0], index=ion_index[:1])
    property_instance = HydrogenContinuumLTEPopulations(None)

    ion_density, level_density = property_instance.calculate(
        number_density=number_density,
        electron_densities=current_electron_densities,
        levels=level_index,
        hydrogen_continuum_partition_function=partition_function,
        thermal_g_electron=np.array([0.5]),
        beta_electron=np.array([1.0]),
        ionization_data=ionization_data,
        thermal_lte_level_boltzmann_factor=thermal_level_factor,
        thermal_lte_partition_function=partition_function,
    )

    pdt.assert_frame_equal(
        ion_density,
        pd.DataFrame([[10.0], [2.0]], index=ion_index, columns=[0]),
    )
    pdt.assert_frame_equal(
        level_density,
        pd.DataFrame([[20.0], [2.0]], index=level_index),
    )


def test_iip_factory_configures_public_stimulated_emission_factor() -> None:
    class AtomData:
        photoionization_data = object()
        selected_atomic_numbers = pd.Index([1, 6])

        def prepare_atom_data(self, *_args: object, **_kwargs: object) -> None:
            return None

    atom_data = AtomData()
    factory = IIPPlasmaSolverFactory(atom_data)
    factory.continuum_interaction_species = []
    nlte_species = [(1, 0), (6, 1)]

    factory.prepare_factory(
        [1, 6],
        "tardis.plasma.properties.iip_property_collections",
        nlte_species,
    )

    assert factory.property_kwargs[StimulatedEmissionFactor] == {
        "nlte_species": set(nlte_species)
    }


def test_coupled_state_property_publishes_nebular_state() -> None:
    columns = pd.Index([0], name="shell")
    ion_index = pd.MultiIndex.from_tuples(
        [(6, ion_number) for ion_number in range(7)],
        names=["atomic_number", "ion_number"],
    )
    level_index = pd.MultiIndex.from_tuples(
        [], names=["atomic_number", "ion_number", "level_number"]
    )
    ionization_index = pd.MultiIndex.from_tuples(
        [(6, ion_number) for ion_number in range(6)],
        names=["atomic_number", "ion_number"],
    )
    phi_index = pd.MultiIndex.from_tuples(
        [(6, ion_number) for ion_number in range(1, 7)],
        names=["atomic_number", "ion_number"],
    )

    previous_beta = pd.DataFrame(columns=columns)
    previous_electrons = pd.Series([7.0], index=columns)
    previous_levels = pd.DataFrame(index=level_index, columns=columns)
    previous_ions = pd.DataFrame(
        np.ones((7, 1)), index=ion_index, columns=columns
    )
    previous_copies = (
        previous_beta.copy(deep=True),
        previous_electrons.copy(deep=True),
        previous_levels.copy(deep=True),
        previous_ions.copy(deep=True),
    )
    number_density = pd.DataFrame(
        [[1.0]],
        index=pd.Index([6], name="atomic_number"),
        columns=columns,
    )
    level_factor = pd.DataFrame(index=level_index, columns=columns)
    thermal_partition_function = pd.DataFrame(
        1.0, index=ion_index, columns=columns
    )
    phi = pd.DataFrame(3.0, index=phi_index, columns=columns)
    property_instance = IIPContinuumPopulations(None)

    result = property_instance.calculate(
        atomic_data=None,
        nlte_data=SimpleNamespace(nlte_species=set()),
        continuum_interaction_species=pd.MultiIndex.from_tuples(
            [], names=["atomic_number", "ion_number"]
        ),
        t_electrons=np.array([10_000.0]),
        dilute_planckian_radiation_field=object(),
        j_blues=None,
        previous_beta_sobolev=previous_beta,
        previous_electron_densities=previous_electrons,
        previous_level_number_density=previous_levels,
        previous_ion_number_density=previous_ions,
        number_density=number_density,
        general_level_boltzmann_factor=level_factor,
        thermal_g_electron=np.array([0.5]),
        beta_electron=np.array([1.0]),
        thermal_lte_level_boltzmann_factor=level_factor,
        thermal_lte_partition_function=thermal_partition_function,
        ionization_data=pd.Series(0.0, index=ionization_index),
        phi=phi,
        g=pd.Series(dtype=float, index=level_index),
        lines_lower_level_index=np.array([], dtype=int),
        lines_upper_level_index=np.array([], dtype=int),
        metastability=pd.Series(dtype=bool, index=level_index),
        lines=pd.DataFrame(),
        time_explosion=1 * u.day,
        photoionization_rate_estimator=None,
        stimulated_recombination_rate_estimator=None,
    )

    pdt.assert_frame_equal(result[0], level_factor)
    pdt.assert_frame_equal(
        result[1],
        pd.DataFrame(
            index=pd.MultiIndex.from_tuples(
                [], names=["atomic_number", "ion_number"]
            ),
            columns=columns,
        ),
    )
    pdt.assert_frame_equal(
        result[2],
        pd.DataFrame(
            np.full((7, 1), 1.0 / 7.0), index=ion_index, columns=columns
        ),
    )
    pdt.assert_series_equal(result[3], pd.Series([3.0], index=columns))
    pdt.assert_frame_equal(result[4], level_factor)
    pdt.assert_frame_equal(previous_beta, previous_copies[0])
    pdt.assert_series_equal(previous_electrons, previous_copies[1])
    pdt.assert_frame_equal(previous_levels, previous_copies[2])
    pdt.assert_frame_equal(previous_ions, previous_copies[3])


def test_estimator_update_requires_complete_previous_population_state() -> None:
    """Estimator updates must not mix current rates with partial old populations."""
    level_index = pd.MultiIndex.from_tuples(
        [], names=["atomic_number", "ion_number", "level_number"]
    )
    ion_index = pd.MultiIndex.from_tuples(
        [], names=["atomic_number", "ion_number"]
    )

    with pytest.raises(ValueError, match="complete previous population state"):
        IIPContinuumPopulations(None).calculate(
            atomic_data=None,
            nlte_data=None,
            continuum_interaction_species=ion_index,
            t_electrons=np.array([]),
            dilute_planckian_radiation_field=None,
            j_blues=pd.DataFrame(),
            previous_beta_sobolev=pd.DataFrame(),
            previous_electron_densities=None,
            previous_level_number_density=pd.DataFrame(index=level_index),
            previous_ion_number_density=pd.DataFrame(index=ion_index),
            number_density=pd.DataFrame(),
            general_level_boltzmann_factor=pd.DataFrame(index=level_index),
            thermal_g_electron=np.array([]),
            beta_electron=np.array([]),
            thermal_lte_level_boltzmann_factor=pd.DataFrame(index=level_index),
            thermal_lte_partition_function=pd.DataFrame(index=ion_index),
            ionization_data=pd.Series(dtype=float),
            phi=pd.DataFrame(),
            g=pd.Series(dtype=float),
            lines_lower_level_index=np.array([], dtype=int),
            lines_upper_level_index=np.array([], dtype=int),
            metastability=pd.Series(dtype=bool),
            lines=pd.DataFrame(),
            time_explosion=1 * u.day,
            photoionization_rate_estimator=object(),
            stimulated_recombination_rate_estimator=object(),
        )
