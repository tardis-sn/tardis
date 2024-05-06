from copy import deepcopy

import numpy.testing as npt
import pandas.testing as pdt
import pytest

from tardis.transport.montecarlo.estimators.continuum_radfield_properties import (
    DiluteBlackBodyContinuumPropertiesSolver,
    MCContinuumPropertiesSolver,
)
from tardis.simulation import Simulation


def test_continuum_estimators(
    continuum_config,
    nlte_atomic_dataset,
):
    nlte_atomic_dataset = deepcopy(nlte_atomic_dataset)
    continuum_simulation = Simulation.from_config(
        continuum_config,
        atom_data=nlte_atomic_dataset,
        virtual_packet_logging=False,
    )
    # continuum_simulation.run_convergence()
    continuum_properties_solver_dilute_bb = (
        DiluteBlackBodyContinuumPropertiesSolver(
            continuum_simulation.plasma.atomic_data
        )
    )

    continuum_properties_dilute_bb = (
        continuum_properties_solver_dilute_bb.solve(
            continuum_simulation.simulation_state.radiation_field_state,
            continuum_simulation.plasma.t_electrons,
        )
    )

    continuum_plasma = continuum_simulation.plasma

    pdt.assert_frame_equal(
        continuum_properties_dilute_bb.photo_ionization_rate_coefficient,
        continuum_simulation.plasma.gamma,
    )
    stimulated_recomb_rate_coeff = (
        continuum_properties_dilute_bb.stimulated_recomb_rate_factor
        * continuum_plasma.phi_ik.loc[
            continuum_properties_dilute_bb.stimulated_recomb_rate_factor.index
        ]
    )
    pdt.assert_frame_equal(
        stimulated_recomb_rate_coeff, continuum_plasma.alpha_stim
    )

    continuum_simulation.run_final()

    continuum_properties_solver_mc = MCContinuumPropertiesSolver(
        continuum_simulation.plasma.atomic_data
    )
    transport_state = continuum_simulation.transport.transport_state
    continuum_properties_mc = continuum_properties_solver_mc.solve(
        transport_state.radfield_mc_estimators,
        transport_state.time_of_simulation,
        transport_state.geometry_state.volume,
    )

    continuum_plasma.update(
        gamma_estimator=transport_state.radfield_mc_estimators.photo_ion_estimator,
        alpha_stim_estimator=transport_state.radfield_mc_estimators.stim_recomb_estimator,
        bf_heating_coeff_estimator=transport_state.radfield_mc_estimators.bf_heating_estimator,
        stim_recomb_cooling_coeff_estimator=transport_state.radfield_mc_estimators.stim_recomb_cooling_estimator,
    )

    pdt.assert_frame_equal(
        continuum_properties_mc.photo_ionization_rate_coefficient,
        continuum_simulation.plasma.gamma,
    )
    stimulated_recomb_rate_coeff = (
        continuum_properties_mc.stimulated_recomb_rate_factor
        * continuum_plasma.phi_ik.loc[
            continuum_properties_dilute_bb.stimulated_recomb_rate_factor.index
        ]
    )
    pdt.assert_frame_equal(
        stimulated_recomb_rate_coeff, continuum_plasma.alpha_stim
    )
