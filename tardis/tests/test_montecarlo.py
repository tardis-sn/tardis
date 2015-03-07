import numpy as np
import yaml
import pytest
import os
from numpy import testing
from tardis import montecarlo
from tardis.base import run_tardis
from tardis.io.config_reader import Configuration

test_line_list = np.array([10, 9, 8, 7, 6, 5, 5, 4, 3, 2, 1]).astype(np.float64)

@pytest.mark.skipif(not pytest.config.getvalue("atomic-dataset"),
                    reason='--atomic_database was not specified')

class TestMonteCarlo():

    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(self):
        self.atom_data_filename = os.path.expanduser(os.path.expandvars(
            pytest.config.getvalue('atomic-dataset')))
        assert os.path.exists(self.atom_data_filename), ("{0} atomic datafiles "
                                                         "does not seem to "
                                                         "exist".format(
            self.atom_data_filename))
        self.config_yaml = yaml.load(open('tardis/io/tests/data/tardis_configv1_verysimple.yml'))
        self.config_yaml['atom_data'] = self.atom_data_filename

        self.model = run_tardis(self.config_yaml)

    def test_montecarlo_radial1d():
        output_nu, output_energies, j_estimators, nubar_estimators, \
        last_line_interaction_in_id, last_line_interaction_out_id, \
        last_interaction_type, last_line_interaction_shell_id = \
            montecarlo.montecarlo.montecarlo_radial1d(self.model, 0)
#        testing.assert_almost_equal(ouput_nu[0], 494994844798236.5)
#        testing.assert_almost_equal(output_energies[0], 2.0570649803508353e-06)
#        testing.assert_almost_equal(j_estimators[0], 114564286557777.5)
#        testing.assert_almost_equal(nubar_estimators[0], 1.0088596598917336e+29)
#        assert last_line_interaction_in_id[0] == -1
#        assert last_line_interaction_out_id[0] == -1
#        assert last_interaction_type[0] == -1
#        assert last_line_interaction_shell_id[0] == -1

# @pytest.mark.parametrize(("insert_value", "expected_insert_position"), [
#     (9.5, 0),
#     (8.5, 1),
#     (7.5, 2),
#     (6.5, 3),
#     (5.5, 4),
#     (5.2, 4),
#     (4.5, 6),
#     (3.5, 7),
#     (2.5, 8),
#     (1.5, 9)])
# def test_binary_search(insert_value, expected_insert_position):
#     insert_position = montecarlo.binary_search_wrapper(test_line_list, insert_value, 0, len(test_line_list))
#     assert insert_position == expected_insert_position


# @pytest.mark.parametrize(("insert_value"), [
#     (10.5),
#     (0.5)])
# def test_binary_search_out_of_bounds(insert_value, capsys):
#     with pytest.raises(ValueError):
#         insert_position = montecarlo.binary_search_wrapper(test_line_list, insert_value, 0, len(test_line_list)-1)

# @pytest.mark.parametrize(("insert_value", "expected_insert_position"), [
#     (10.5, 0),
#     (0.5, len(test_line_list))])
# def test_line_search_out_of_bounds(insert_value, expected_insert_position):
#     insert_position = montecarlo.line_search_wrapper(test_line_list,
#                             insert_value, len(test_line_list))
#     assert insert_position == expected_insert_position

# def test_compute_distance2outer():
#     assert montecarlo.compute_distance2outer_wrapper(0.0, 0.5, 1.0) == 1.0
#     assert montecarlo.compute_distance2outer_wrapper(1.0, 0.5, 1.0) == 0.0
#     assert montecarlo.compute_distance2outer_wrapper(0.3, 1.0, 1.0) == 0.7
#     assert montecarlo.compute_distance2outer_wrapper(0.3, -1.0, 1.0) == 1.3
#     assert montecarlo.compute_distance2outer_wrapper(0.5, 0.0, 1.0) == np.sqrt(0.75)

# def test_compute_distance2inner():
#     assert montecarlo.compute_distance2inner_wrapper(1.5, -1.0, 1.0) == 0.5
#     assert montecarlo.compute_distance2inner_wrapper(0.0, 0.0, 0.0) == montecarlo.miss_distance
#     assert montecarlo.compute_distance2inner_wrapper(1.2, -0.7, 1.0) == 0.3246360509309949

# def test_compute_distance2line():
#     assert montecarlo.compute_distance2line_wrapper(2.20866912e+15, -0.251699059004, 1.05581082105e+15, 1.06020910733e+15, 1693440.0, 5.90513983371e-07, 1.0602263591e+15, 1.06011723237e+15, 2) == 344430881691490.5
#     assert montecarlo.compute_distance2line_wrapper(2.23434667994e+15, -0.291130548401, 1.05581082105e+15, 1.06733618121e+15, 1693440.0, 5.90513983371e-07, 1.06738407486e+15, 1.06732933961e+15, 3) == 96296282395637.2
#     with pytest.raises(RuntimeError):
#         montecarlo.compute_distance2line_wrapper(1.0, 1.0, 1.0, 10.0, 15.0, 1.0 / 15.0, 0.0, 0.0, 0)

# def test_compute_distance2electron():
#     assert montecarlo.compute_distance2electron_wrapper(0.0, 0.0, 2.0, 2.0) == 4.0

