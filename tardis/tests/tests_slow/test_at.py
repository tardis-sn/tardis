import pytest
from base import SlowTest


@pytest.mark.skipif(not pytest.config.getoption("--run-slow"),
                    reason='this is a slow test, add --run-slow to run')
class TestAT(SlowTest):
    """
    Integration test for a run with the AT setup.
    """

    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(self):
        """
        This method does an initial setup of making the configuration and performing
        a single run of this integration test from AT setup.
        """

        super(TestAT, self).base_setup(
            model_dir='abn_tom_test',
            config_file='abn_tom_test.yml',
            abundances='abn_tom_test_abundances.dat',
            densities='abn_tom_test_densities.dat')

    def test_j_estimators(self):
        super(TestAT, self).test_j_estimators()

    def test_j_blue_estimators(self):
        super(TestAT, self).test_j_blue_estimators()

    def test_last_line_interactions(self):
        super(TestAT, self).test_last_line_interactions()

    def test_nubar_estimators(self):
        super(TestAT, self).test_nubar_estimators()

    def test_ws(self):
        super(TestAT, self).test_ws()

    def test_spectrum(self):
        super(TestAT, self).test_spectrum()

    def test_montecarlo_properties(self):
        super(TestAT, self).test_montecarlo_properties()

    def test_shell_temperature(self):
        super(TestAT, self).test_shell_temperature()
