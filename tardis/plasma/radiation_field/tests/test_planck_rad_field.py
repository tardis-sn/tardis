from types import SimpleNamespace

import numpy as np
from astropy import units as u
from numpy.testing import assert_array_almost_equal

from tardis.model import SimulationState
from tardis.plasma.radiation_field import DilutePlanckianRadiationField


def test_dilution_factor_setter_accepts_read_only_initial_array() -> None:
    geometry = SimpleNamespace(
        v_inner_boundary_index=1,
        v_outer_boundary_index=3,
        no_of_shells_active=2,
    )
    dilution_factor = np.array([0.1, 0.2, 0.3, 0.4])
    dilution_factor.flags.writeable = False
    radiation_field_state = DilutePlanckianRadiationField(
        np.ones(4) * 10000 * u.K,
        dilution_factor,
        geometry,
    )
    simulation_state = SimulationState(
        geometry=geometry,
        composition=None,
        radiation_field_state=radiation_field_state,
        time_explosion=None,
        packet_source=None,
    )

    simulation_state.dilution_factor = np.array([0.5, 0.6])

    assert_array_almost_equal(
        simulation_state.radiation_field_state.dilution_factor,
        np.array([0.1, 0.5, 0.6, 0.4]),
    )
