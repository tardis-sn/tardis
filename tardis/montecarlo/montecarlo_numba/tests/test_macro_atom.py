import pytest
import tardis.montecarlo.montecarlo_numba.macro_atom as macro_atom
import numpy as np
from numba import njit

@pytest.mark.parametrize(
    'expected',
    [5259]
    )
def test_macro_atom(static_packet, verysimple_numba_plasma, verysimple_numba_model, expected):
    static_packet.initialize_line_id(verysimple_numba_plasma, verysimple_numba_model)
    result = macro_atom.macro_atom(static_packet, verysimple_numba_plasma)
    assert result == expected