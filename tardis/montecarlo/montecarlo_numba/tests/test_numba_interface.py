import pytest
import tardis.montecarlo.montecarlo_numba.numba_interface as numba_interface
import numpy.testing as npt
import numpy as np

@pytest.mark.parametrize( 
    'input_params', 
    ['scatter', 'macroatom', 'downbranch']
)
def test_numba_plasma_initialize(nb_simulation_verysimple, input_params):
    line_interaction_type = input_params
    plasma = nb_simulation_verysimple.plasma
    actual = numba_interface.numba_plasma_initialize(plasma, line_interaction_type)
    
    npt.assert_allclose(actual.electron_density, plasma.electron_densities.values)
    npt.assert_allclose(actual.line_list_nu, plasma.atomic_data.lines.nu.values)
    npt.assert_allclose(actual.tau_sobolev, plasma.tau_sobolevs.values)
    if line_interaction_type == 'scatter':
        empty = np.zeros(1, dtype=np.int64)
        npt.assert_allclose(actual.transition_probabilities, np.zeros((1, 1), dtype=np.float64))
        npt.assert_allclose(actual.line2macro_level_upper, empty)
        npt.assert_allclose(actual.macro_block_references, empty)
        npt.assert_allclose(actual.transition_type, empty)
        npt.assert_allclose(actual.destination_level_id, empty)
        npt.assert_allclose(actual.transition_line_id, empty)
    else:
        npt.assert_allclose(actual.transition_probabilities, plasma.transition_probabilities.values)
        npt.assert_allclose(actual.line2macro_level_upper, plasma.atomic_data.lines_upper2macro_reference_idx)
        npt.assert_allclose(actual.macro_block_references, plasma.atomic_data.macro_atom_references[
            'block_references'].values)
        npt.assert_allclose(actual.transition_type, plasma.atomic_data.macro_atom_data[
            'transition_type'].values)
        npt.assert_allclose(actual.destination_level_id, plasma.atomic_data.macro_atom_data[
            'destination_level_idx'].values)
        npt.assert_allclose(actual.transition_line_id, plasma.atomic_data.macro_atom_data[
            'lines_idx'].values)


@pytest.mark.xfail(reason='To be implemented')
def test_configuration_initialize():
    assert False