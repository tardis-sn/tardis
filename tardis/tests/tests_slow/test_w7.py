import pytest
from base import SlowTest


@pytest.mark.skipif(not pytest.config.getoption("--run-slow"),
                    reason='this is a slow test, add --run-slow to run')
class TestW7(SlowTest):
    """
    Integration test for a run with the Stratified W7 setup.
    """

    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(self):
        """
        This method does an initial setup of making the configuration and performing
        a single run of this integration test from Stratified W7 setup.
        """

        super(TestW7, self).base_setup(
            model_dir='w7_13d',
            config_file='tardis_w7.yml',
            abundances='tardis_w7_13d_abundances.dat',
            densities='tardis_w7_13d_densities.dat')

    def test_j_estimators(self):
        super(TestW7, self).test_j_estimators()

    def test_j_blue_estimators(self):
        super(TestW7, self).test_j_blue_estimators()

    def test_last_line_interactions(self):
        super(TestW7, self).test_last_line_interactions()

    def test_nubar_estimators(self):
        super(TestW7, self).test_nubar_estimators()

    def test_ws(self):
        super(TestW7, self).test_ws()

    def test_spectrum(self):
        super(TestW7, self).test_ws()

    def test_montecarlo_properties(self):
        super(TestW7, self).test_montecarlo_properties()

    def test_shell_temperature(self):
        super(TestW7, self).test_shell_temperature()
