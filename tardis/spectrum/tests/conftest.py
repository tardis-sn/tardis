import pytest

from tardis.opacities.opacity_state import opacity_state_initialize


@pytest.fixture(scope="package")
def simulation_verysimple_opacity_state(simulation_verysimple):
    return opacity_state_initialize(
        simulation_verysimple.plasma,
        line_interaction_type="macroatom",
        disable_line_scattering=False,
        continuum_processes_enabled=False,
    )
