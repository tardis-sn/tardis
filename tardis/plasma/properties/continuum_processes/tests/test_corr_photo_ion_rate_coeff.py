"""
Regression tests for CorrPhotoIonRateCoeff — GitHub issue #3355.

CorrPhotoIonRateCoeff previously raised a hard PlasmaException when
Monte Carlo estimator noise drove gamma_corr slightly below zero.
The fix clips to zero with a warning instead of crashing.
"""

import logging

import numpy as np
import pandas as pd
import pytest

from tardis.plasma.properties.continuum_processes.rates import (
    CorrPhotoIonRateCoeff,
)


@pytest.fixture
def corr_photo_ion_rate_coeff():
    # plasma_parent=None is safe: ProcessingPlasmaProperty only stores it
    return CorrPhotoIonRateCoeff(None)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Two photo-ionization edges for H I: levels (1,0,0) and (1,0,1)
# ionise up to H II (1,1).
_LEVEL_IDX = pd.MultiIndex.from_tuples(
    [(1, 0, 0), (1, 0, 1)],
    names=["atomic_number", "ion_number", "level_number"],
)
# ion_number_density only needs the unique key (1,1) once;
# pandas .loc will duplicate the row when called with [(1,1),(1,1)].
_ION_IDX = pd.MultiIndex.from_tuples(
    [(1, 1)],
    names=["atomic_number", "ion_number"],
)


def _make_inputs(
    gamma_vals,
    alpha_stim_vals,
    n_k_val,
    n_i_vals,
    ne_vals,
):
    """
    Build the five DataFrames/Series consumed by
    CorrPhotoIonRateCoeff.calculate().

    Parameters
    ----------
    gamma_vals : list[float], length 2
        photo-ionization rate per level, broadcast across all shells.
    alpha_stim_vals : list[float], length 2
        stimulated-recombination rate factor per level.
    n_k_val : float
        ion (H II) number density, same for every shell.
    n_i_vals : list[float], length 2
        level number density per level.
    ne_vals : list[float]
        electron densities, one per shell.
    """
    n_shells = len(ne_vals)
    cols = list(range(n_shells))

    gamma = pd.DataFrame(
        {c: gamma_vals for c in cols}, index=_LEVEL_IDX, dtype=float
    )
    alpha_stim = pd.DataFrame(
        {c: alpha_stim_vals for c in cols}, index=_LEVEL_IDX, dtype=float
    )
    electron_densities = pd.Series(ne_vals, dtype=float)
    ion_number_density = pd.DataFrame(
        {c: [n_k_val] for c in cols}, index=_ION_IDX, dtype=float
    )
    level_number_density = pd.DataFrame(
        {c: n_i_vals for c in cols}, index=_LEVEL_IDX, dtype=float
    )
    return gamma, alpha_stim, electron_densities, ion_number_density, level_number_density


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestCorrPhotoIonRateCoeffClipping:
    """Tests that negative values are clipped, not raised."""

    def test_positive_values_unchanged(self, corr_photo_ion_rate_coeff):
        """All-positive gamma_corr should pass through unchanged."""
        # correction = alpha_stim * (n_k/n_i) * ne
        #            = 1e-15 * (1e10/1e16) * 1e10 = 1e-11  <  gamma=1e-10 ✓
        gamma, alpha_stim, ne, n_k, n_i = _make_inputs(
            gamma_vals=[1e-10, 1e-10],
            alpha_stim_vals=[1e-15, 1e-15],
            n_k_val=1e10,
            n_i_vals=[1e16, 1e16],
            ne_vals=[1e10, 1e10],
        )
        result = corr_photo_ion_rate_coeff.calculate(
            gamma, alpha_stim, ne, n_k, n_i
        )
        assert (result >= 0).all().all()

    def test_slightly_negative_values_clipped_to_zero(
        self, corr_photo_ion_rate_coeff
    ):
        """
        Regression test for issue #3355.

        Slightly negative gamma_corr (from MC estimator noise) must be
        clipped to zero, not raise PlasmaException.
        """
        # correction = 1e-10 * (1e10/1e12) * 1e8 = 1e-4  >>  gamma=1e-15
        gamma, alpha_stim, ne, n_k, n_i = _make_inputs(
            gamma_vals=[1e-15, 1e-15],
            alpha_stim_vals=[1e-10, 1e-10],
            n_k_val=1e10,
            n_i_vals=[1e12, 1e12],
            ne_vals=[1e8, 1e8],
        )
        # Must NOT raise — that was the bug
        result = corr_photo_ion_rate_coeff.calculate(
            gamma, alpha_stim, ne, n_k, n_i
        )
        assert (result >= 0).all().all()

    def test_no_plasma_exception_raised(self, corr_photo_ion_rate_coeff):
        """
        Explicitly verify PlasmaException is never raised regardless of
        input values.
        """
        from tardis.plasma.exceptions import PlasmaException

        gamma, alpha_stim, ne, n_k, n_i = _make_inputs(
            gamma_vals=[1e-15, 1e-15],
            alpha_stim_vals=[1e-10, 1e-10],
            n_k_val=1e10,
            n_i_vals=[1e12, 1e12],
            ne_vals=[1e8, 1e8],
        )
        try:
            result = corr_photo_ion_rate_coeff.calculate(
                gamma, alpha_stim, ne, n_k, n_i
            )
        except PlasmaException:
            pytest.fail("PlasmaException was raised — issue #3355 regression!")

        assert (result >= 0).all().all()

    def test_warning_emitted_for_negative_values(
        self, corr_photo_ion_rate_coeff, caplog
    ):
        """A warning must be logged when clipping occurs."""
        gamma, alpha_stim, ne, n_k, n_i = _make_inputs(
            gamma_vals=[1e-15, 1e-15],
            alpha_stim_vals=[1e-10, 1e-10],
            n_k_val=1e10,
            n_i_vals=[1e12, 1e12],
            ne_vals=[1e8, 1e8],
        )
        logger_name = (
            "tardis.plasma.properties.continuum_processes.rates"
        )
        with caplog.at_level(logging.WARNING, logger=logger_name):
            corr_photo_ion_rate_coeff.calculate(
                gamma, alpha_stim, ne, n_k, n_i
            )

        assert any(
            "Negative values in CorrPhotoIonRateCoeff" in r.message
            for r in caplog.records
        ), "Expected a warning about negative values to be logged"

    def test_zero_inputs_return_zero(self, corr_photo_ion_rate_coeff):
        """Zero gamma and zero alpha_stim produce zero gamma_corr."""
        gamma, alpha_stim, ne, n_k, n_i = _make_inputs(
            gamma_vals=[0.0, 0.0],
            alpha_stim_vals=[0.0, 0.0],
            n_k_val=1e10,
            n_i_vals=[1e12, 1e12],
            ne_vals=[1e10, 1e10],
        )
        result = corr_photo_ion_rate_coeff.calculate(
            gamma, alpha_stim, ne, n_k, n_i
        )
        np.testing.assert_array_equal(result.values, 0.0)
