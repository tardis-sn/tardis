"""
Basic TARDIS Benchmark.
"""
import numpy as np
import numpy.testing as npt
import pandas as pd
from asv_runner.benchmarks.mark import parameterize, skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis.plasma.properties.nlte_excitation_data import NLTEExcitationData
from tardis.plasma.properties.nlte_rate_equation_solver import (
    prepare_r_uls_r_lus,
    prepare_bound_bound_rate_matrix,
    create_coll_exc_deexc_matrix,
)


# @skip_benchmark
class BenchmarkPlasmaNlteExcitation(BenchmarkBase):
    """
    Class to benchmark the NLTE excitation function.
    """

    def time_prepare_bound_bound_rate_matrix(self):
        nlte_atomic_dataset = self.nlte_atomic_dataset
        # TODO: Needs to work in the class RegressionData in BenchmarkBase
        custom_request = self.CustomPyTestRequest(
            tardis_regression_data_path="/app/tardis-regression-data",
            node_name="test_prepare_bound_bound_rate_matrix",
            node_module_name="tardis.plasma.tests.test_nlte_excitation",
            regression_data_dir="tardis/plasma/tests/test_nlte_excitation",
        )
        regression_data = self.regression_data(custom_request)

        simple_excitation_species = [(1, 0)]
        copy_atomic_dataset = nlte_atomic_dataset
        copy_atomic_dataset.levels = nlte_atomic_dataset.levels[
            nlte_atomic_dataset.levels.index.get_level_values("level_number") < 5
            ]
        lines_filtered = nlte_atomic_dataset.lines[
            nlte_atomic_dataset.lines.index.get_level_values("level_number_lower")
            < 5
            ]
        copy_atomic_dataset.lines = lines_filtered[
            lines_filtered.index.get_level_values("level_number_upper") < 5
            ]
        simple_nlte_data = NLTEExcitationData(
            copy_atomic_dataset.lines, simple_excitation_species
        )
        simple_number_of_shells = 1
        simple_j_blues_matrix = np.linspace(
            0.1, 1.0, copy_atomic_dataset.lines.index.size
        )
        simple_j_blues = pd.DataFrame(
            simple_j_blues_matrix,
            index=copy_atomic_dataset.lines.index,
            columns=["0"],
        )
        simple_number_of_levels = copy_atomic_dataset.levels.energy.loc[
            simple_excitation_species[0]
        ].count()
        (
            lines_index,
            r_ul_index,
            r_ul_matrix,
            r_lu_index,
            r_lu_matrix,
        ) = prepare_r_uls_r_lus(
            simple_number_of_levels,
            simple_number_of_shells,
            simple_j_blues,
            simple_excitation_species[0],
            simple_nlte_data,
        )
        simple_beta_sobolev_matrix = np.linspace(
            2.5, 4.7, copy_atomic_dataset.lines.index.size
        )
        simple_beta_sobolev = pd.DataFrame(
            simple_beta_sobolev_matrix,
            index=copy_atomic_dataset.lines.index,
            columns=["0"],
        )
        actual_rate_matrix = prepare_bound_bound_rate_matrix(
            simple_number_of_levels,
            lines_index,
            r_ul_index,
            r_ul_matrix,
            r_lu_index,
            r_lu_matrix,
            simple_beta_sobolev,
        )
        # if this test fails the first thing to check is if the reshape in the
        # methods made a view or a copy. If it's a copy rewrite the function.
        # TODO: allow rtol=1e-6
        expected_rate_matrix = regression_data.sync_ndarray(actual_rate_matrix)
        npt.assert_allclose(actual_rate_matrix, expected_rate_matrix, rtol=1e-6)

    @parameterize({'Matrix': [
        {
            "coll_exc_coeff_values": [1, -2, 3],
            "coll_deexc_coeff_values": [4, 9, 10],
            "number_of_levels": 3,
        },
        {
            "coll_exc_coeff_values": [0.21, 0.045, 0.1234],
            "coll_deexc_coeff_values": [0.7865, 0.987, 0.00123],
            "number_of_levels": 3,
        },
    ]})
    def time_coll_exc_deexc_matrix(self, matrix):
        """
        Checks the NLTERateEquationSolver.create_coll_exc_deexc_matrix for simple values of species with 3 levels.
        NOTE: Values used for testing are not physical.
        """
        coll_exc_coeff_values = matrix["coll_exc_coeff_values"]
        coll_deexc_coeff_values = matrix["coll_deexc_coeff_values"]
        number_of_levels = matrix["number_of_levels"]
        index = pd.MultiIndex.from_tuples(
            [(0, 1), (0, 2), (1, 2)],
            names=["level_number_lower", "level_number_upper"],
        )
        exc_coeff = pd.Series(coll_exc_coeff_values, index=index)
        deexc_coeff = pd.Series(coll_deexc_coeff_values, index=index)
        obtained_coeff_matrix = create_coll_exc_deexc_matrix(
            exc_coeff, deexc_coeff, number_of_levels
        )
        # TODO: Needs to work in the class RegressionData in BenchmarkBase
        test_parameter_id = 0
        if coll_exc_coeff_values != [1, -2, 3]:
            test_parameter_id = 1
        custom_request = self.CustomPyTestRequest(
            "/app/tardis-regression-data",
            f"test_coll_exc_deexc_matrix__coll_exc_coeff_values{test_parameter_id}-coll_deexc_coeff_values{test_parameter_id}-3__",
            "tardis.plasma.tests.test_nlte_excitation",
            "tardis/plasma/tests/test_nlte_excitation",
        )
        regression_data = self.regression_data(custom_request)

        expected_obtained_coeff_matrix = regression_data.sync_ndarray(
            obtained_coeff_matrix
        )
        npt.assert_allclose(expected_obtained_coeff_matrix, obtained_coeff_matrix)
