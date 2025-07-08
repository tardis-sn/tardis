import pytest

from tardis.spectrum.formal_integral.base import check
from tardis.transport.montecarlo.configuration import montecarlo_globals

@pytest.mark.parametrize(
        "line_interaction_type", 
            ("downbranch", 
             "macroatom", 
            pytest.param("?", marks=pytest.mark.xfail)
            )

)
def test_check(simulation_verysimple, line_interaction_type):

    sim_state = simulation_verysimple.simulation_state
    plasma = simulation_verysimple.plasma
    transport = simulation_verysimple.transport
    transport.line_interaction_type = line_interaction_type

    assert check(sim_state, plasma, transport)
)
    assert not check(None, plasma, transport, raises=False)



