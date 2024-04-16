import os

from astropy import units as u
import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

from tardis.transport.montecarlo.packet_source import (
    BlackBodySimpleSource,
    BlackBodySimpleSourceRelativistic,
)
from tardis.transport.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)
from tardis.tests.fixtures.regression_data import RegressionData


@pytest.fixture(scope="function")
def packet_unit_test_fpath(tardis_ref_path):
    """
    Path to `packet_unittest.h5`.

    Parameters
    ----------
    tardis_ref_path : pd.HDFStore

    Returns
    -------
    os.path
    """
    return os.path.abspath(
        os.path.join(tardis_ref_path, "packet_unittest.h5")
    )

class TestBlackBodySimpleSource:
    @pytest.fixture(scope="function")
    def blackbodysimplesource(self, request, regression_data):
        """
        Create BlackBodySimpleSource instance.

        Yields
        -------
        tardis.transport.montecarlo.packet_source.BlackBodySimpleSource
        """
        cls = type(self)
        montecarlo_configuration.LEGACY_MODE_ENABLED = True
        bb = BlackBodySimpleSource(
            radius=123,
            temperature = 10000 * u.K,
            base_seed=1963, 
            legacy_second_seed=2508
        )
        nus = bb.create_packet_nus(100).value
        mus = bb.create_packet_mus(100)
        unif_energies = bb.create_packet_energies(100).value
        regression_data.sync_hdf_store(bb)

        cls.regression_data =  regression_data
        cls.bb = bb

        yield nus, mus, unif_energies
        montecarlo_configuration.LEGACY_MODE_ENABLED = False

    def test_bb_packet_sampling(
        self,
        request,
        tardis_ref_data,
        packet_unit_test_fpath,
        blackbodysimplesource,
    ):
        """
        Parameters
        ----------
        request : _pytest.fixtures.SubRequest
        tardis_ref_data: pd.HDFStore
        packet_unit_test_fpath: os.path
        """
        if request.config.getoption("--generate-reference"):
            ref_bb = pd.read_hdf(packet_unit_test_fpath, key="/blackbody")
            ref_bb.to_hdf(
                tardis_ref_data, key="/packet_unittest/blackbody", mode="a"
            )
            pytest.skip("Reference data was generated during this run.")

        ref_df = tardis_ref_data["/packet_unittest/blackbody"]
        nus, mus, unif_energies = blackbodysimplesource
        assert np.all(np.isclose(nus, ref_df["nus"]))
        assert np.all(np.isclose(mus, ref_df["mus"]))
        assert np.all(np.isclose(unif_energies, ref_df["energies"]))

        expected_bb = pd.read_hdf(
            self.regression_data.fpath, key="/black_body_simple_source/scalars"
        )
        assert_allclose(expected_bb.base_seed, self.bb.base_seed)
        assert_allclose(expected_bb.temperature, self.bb.temperature.value)
        assert_allclose(expected_bb.radius, self.bb.radius)
        
class TestBlackBodySimpleSourceRel:
    @pytest.fixture(scope="function")
    def blackbody_simplesource_relativistic(self, request, regression_data):
        """
        Create BlackBodySimpleSourceRelativistic instance.

        Yields
        -------
        tardis.montecarlo.packet_source.BlackBodySimpleSourceRelativistic
        """
        cls = type(self)
        montecarlo_configuration.LEGACY_MODE_ENABLED = True

        bb_rel = BlackBodySimpleSourceRelativistic(
            time_explosion=1123187
        )

        bb_rel.temperature = 10000 * u.K
        bb_rel.beta = 0.25
        nus = bb_rel.create_packet_nus(100).value
        unif_energies = (
            bb_rel.create_packet_energies(
                100
            ).value
        )
        bb_rel._reseed(2508)
        mus = bb_rel.create_packet_mus(10)
        gamma = np.sqrt(1 - bb_rel.beta**2) ** -1

        regression_data.sync_hdf_store(bb_rel)
        cls.regression_data = regression_data
        cls.bb_rel = bb_rel

        yield nus, mus, unif_energies, gamma
        montecarlo_configuration.LEGACY_MODE_ENABLED = False

    def test_bb_packet_sampling_relativistic(
        self,
        request,
        tardis_ref_data,
        blackbody_simplesource_relativistic,
    ):
        """
        Parameters
        ----------
        tardis_ref_data : pd.HDFStore
        blackbody_simplesource_relativistic : tardis.transport.montecarlo.packet_source.BlackBodySimpleSourceRelativistic
        """
        if request.config.getoption("--generate-reference"):
            ref_bb = pd.read_hdf(packet_unit_test_fpath, key="/blackbody")
            ref_bb.to_hdf(
                tardis_ref_data, key="/packet_unittest/blackbody_rel", mode="a"
            )
            pytest.skip("Reference data was generated during this run.")

        nus, mus, unif_energies, gamma = blackbody_simplesource_relativistic
        ref_df = tardis_ref_data["/packet_unittest_rel/blackbody"]
        expected_nus = ref_df["nus"]
        expected_unif_energies = ref_df["energies"] * 1.6 / gamma
        expected_mus = np.array(
            [
                0.60420546,
                0.49899691,
                0.69583288,
                0.96812652,
                0.01544154,
                0.93562304,
                0.44306545,
                0.77010037,
                0.896973,
                0.67876489,
            ]
        )
        expected_bb = pd.read_hdf(
            self.regression_data.fpath, key="/black_body_simple_source/scalars"
        )
        assert_allclose(nus, expected_nus)
        assert_allclose(unif_energies, expected_unif_energies)
        assert_allclose(mus, expected_mus, rtol=1e-6)

        assert_allclose(
            expected_bb.time_explosion, self.bb_rel.time_explosion
        )
        assert_allclose(
            expected_bb.temperature, self.bb_rel.temperature.value
        )

