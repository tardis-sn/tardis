import os

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
        montecarlo_configuration.LEGACY_MODE_ENABLED = True
        cls.bb = BlackBodySimpleSource(base_seed=1963, legacy_second_seed=2508)
        yield cls.bb
        montecarlo_configuration.LEGACY_MODE_ENABLED = False

    @pytest.fixture(scope="class")
    def blackbodysimplesourcerelativistic(self, request):
        """
        Create BlackBodySimpleSourceRelativistic instance.

        Yields
        -------
        tardis.montecarlo.packet_source.BlackBodySimpleSourceRelativistic
        """
        montecarlo_configuration.LEGACY_MODE_ENABLED = True
        bb_rel = BlackBodySimpleSourceRelativistic(
            base_seed=1963, legacy_second_seed=2508
        )
        yield bb_rel
        montecarlo_configuration.LEGACY_MODE_ENABLED = False

    def test_bb_packet_sampling(
        self,
        request,
        tardis_ref_data,
        packet_unit_test_fpath,
        blackbodysimplesource,
    ):
        """
        Test generate_plot_mpl method.

        Parameters
        ----------
        request : _pytest.fixtures.SubRequest
        tardis_ref_data: pd.HDFStore
        packet_unit_test_fpath: os.path
        blackbodysimplesource: tardis.montecarlo.packet_source.BlackBodySimpleSource
        """
        if request.config.getoption("--generate-reference"):
            ref_bb = pd.read_hdf(packet_unit_test_fpath, key="/blackbody")
            ref_bb.to_hdf(
                tardis_ref_data, key="/packet_unittest/blackbody", mode="a"
            )
            pytest.skip("Reference data was generated during this run.")

        ref_df = tardis_ref_data["/packet_unittest/blackbody"]
        self.bb.temperature = 10000
        nus = self.bb.create_packet_nus(100)
        mus = self.bb.create_packet_mus(100)
        unif_energies = self.bb.create_packet_energies(100)
        assert np.all(np.isclose(nus, ref_df["nus"]))
        assert np.all(np.isclose(mus, ref_df["mus"]))
        assert np.all(np.isclose(unif_energies, ref_df["energies"]))

    def test_bb_packet_sampling_relativistic(
        self,
        tardis_ref_data,
        blackbodysimplesourcerelativistic,
    ):
        """
        Parameters
        ----------
        tardis_ref_data : pd.HDFStore
        blackbodysimplesourcerelativistic : tardis.montecarlo.packet_source.BlackBodySimpleSourceRelativistic
        """
        blackbodysimplesourcerelativistic.temperature = 10000
        blackbodysimplesourcerelativistic.beta = 0.25

        nus = blackbodysimplesourcerelativistic.create_packet_nus(100)
        unif_energies = (
            blackbodysimplesourcerelativistic.create_packet_energies(100)
        )
        blackbodysimplesourcerelativistic._reseed(2508)
        mus = blackbodysimplesourcerelativistic.create_packet_mus(10)

        gamma = np.sqrt(1 - 0.25**2) ** -1
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
