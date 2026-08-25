import pytest

from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.simulation import Simulation


@pytest.fixture(scope="package")
def simulation_verysimple_opacity_state(
    simulation_verysimple: Simulation,
) -> OpacityStateNumba:
    """Numba opacity state from the very simple simulation"""
    return simulation_verysimple.opacity_state.to_numba(
        simulation_verysimple.macro_atom_state,
        "macroatom",
    )
