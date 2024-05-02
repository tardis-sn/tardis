import os

from astropy import units as u
import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

from tardis.montecarlo.packet_source import (
    BlackBodySimpleSource,
    BlackBodySimpleSourceRelativistic,
)
from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)


class TestPacketSource:
    @pytest.fixture(scope="class")
    def packet_unit_test_fpath(self, tardis_ref_path):
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

    @pytest.fixture(scope="class")
    def blackbodysimplesource(self, request):
        """
        Create BlackBodySimpleSource instance.

        Yields
        -------
        tardis.montecarlo.packet_source.BlackBodySimpleSource
        """
        cls = type(self)
        cls.bb = BlackBodySimpleSource(
            base_seed=1963, legacy_mode_enabled=True, legacy_second_seed=2508
        )
        yield cls.bb

    @pytest.fixture(scope="class")
    def blackbody_simplesource_relativistic(self, request):
        """
        Create BlackBodySimpleSourceRelativistic instance.

        Yields
        -------
        tardis.montecarlo.packet_source.BlackBodySimpleSourceRelativistic
        """
        bb_rel = BlackBodySimpleSourceRelativistic(
            base_seed=1963, legacy_mode_enabled=True, legacy_second_seed=2508
        )
        yield bb_rel

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
        self.bb.temperature = 10000 * u.K
        nus = self.bb.create_packet_nus(100).value
        mus = self.bb.create_packet_mus(100)
        unif_energies = self.bb.create_packet_energies(100).value
        assert np.all(np.isclose(nus, ref_df["nus"]))
        assert np.all(np.isclose(mus, ref_df["mus"]))
        assert np.all(np.isclose(unif_energies, ref_df["energies"]))

    def test_bb_packet_sampling_relativistic(
        self,
        tardis_ref_data,
        blackbody_simplesource_relativistic,
    ):
        """
        Parameters
        ----------
        tardis_ref_data : pd.HDFStore
        blackbody_simplesource_relativistic : tardis.montecarlo.packet_source.BlackBodySimpleSourceRelativistic
        """
        blackbody_simplesource_relativistic.temperature = 10000 * u.K
        blackbody_simplesource_relativistic.beta = 0.25

        nus = blackbody_simplesource_relativistic.create_packet_nus(100).value
        unif_energies = (
            blackbody_simplesource_relativistic.create_packet_energies(
                100
            ).value
        )
        blackbody_simplesource_relativistic._reseed(2508)
        mus = blackbody_simplesource_relativistic.create_packet_mus(10)

        gamma = np.sqrt(1 - blackbody_simplesource_relativistic.beta**2) ** -1
        ref_df = tardis_ref_data["/packet_unittest/blackbody"]
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

        assert_allclose(nus, expected_nus)
        assert_allclose(unif_energies, expected_unif_energies)
        assert_allclose(mus, expected_mus, rtol=1e-6)
