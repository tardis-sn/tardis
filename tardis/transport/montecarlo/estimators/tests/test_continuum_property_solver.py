from copy import deepcopy

import pandas.testing as pdt
import pytest
from astropy import units as u

from tardis.plasma.equilibrium.rates.photoionization_strengths import (
    AnalyticPhotoionizationCoeffSolver,
)
from tardis.simulation import Simulation
from tardis.transport.montecarlo.estimators.continuum_radfield_properties import (
    ContinuumProperties,
    MCContinuumPropertiesSolver,
)


@pytest.mark.skip("Fix once continuum works again")

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
    continuum_properties_solver_dilute_bb = AnalyticPhotoionizationCoeffSolver(
        continuum_simulation.plasma.atomic_data.photoionization_data
    )

    continuum_properties_dilute_bb = ContinuumProperties(
        *continuum_properties_solver_dilute_bb.solve(
            continuum_simulation.simulation_state.radiation_field_state,
            continuum_simulation.plasma.t_electrons * u.K,
        )
    )

    continuum_plasma = continuum_simulation.plasma

    pdt.assert_frame_equal(
        continuum_properties_dilute_bb.photo_ionization_rate_coefficient,
        continuum_simulation.plasma.gamma,
    )
    stimulated_recomb_rate_coeff = (
        continuum_properties_dilute_bb.stimulated_recombination_rate_factor
        * continuum_plasma.phi_ik.loc[
            continuum_properties_dilute_bb.stimulated_recombination_rate_factor.index
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
        transport_state.estimators_continuum,
        transport_state.time_of_simulation,
        transport_state.geometry_state.volume,
    )

    continuum_plasma.update(
        gamma=continuum_properties_mc.photo_ionization_rate_coefficient,
        alpha_stim_factor=continuum_properties_mc.stimulated_recombination_rate_factor,
        bf_heating_coeff_estimator=None,
        stim_recomb_cooling_coeff_estimator=None,
    )

    pdt.assert_frame_equal(
        continuum_properties_mc.photo_ionization_rate_coefficient,
        continuum_simulation.plasma.gamma,
    )
    stimulated_recomb_rate_coeff = (
        continuum_properties_mc.stimulated_recombination_rate_factor
        * continuum_plasma.phi_ik.loc[
            continuum_properties_dilute_bb.stimulated_recombination_rate_factor.index
        ]
    )
    pdt.assert_frame_equal(
        stimulated_recomb_rate_coeff, continuum_plasma.alpha_stim
    )
