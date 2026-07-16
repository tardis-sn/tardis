import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt

from tardis.opacities.tau_sobolev import (
    calculate_beta_sobolev,
    calculate_sobolev_line_opacity,
)


def test_calculate_sobolev_line_opacity(
    nb_simulation_verysimple, regression_data
):
    legacy_plasma = nb_simulation_verysimple.plasma

    actual = calculate_sobolev_line_opacity(
        legacy_plasma.lines,
        legacy_plasma.level_number_density,
        legacy_plasma.time_explosion,
        legacy_plasma.stimulated_emission_factor,
    )
    expected = regression_data.sync_dataframe(actual)
    pdt.assert_frame_equal(actual, expected)


def test_calculate_beta_sobolevs(nb_simulation_verysimple, regression_data):
    legacy_plasma = nb_simulation_verysimple.plasma

    tau_sobolevs = calculate_sobolev_line_opacity(
        legacy_plasma.lines,
        legacy_plasma.level_number_density,
        legacy_plasma.time_explosion,
        legacy_plasma.stimulated_emission_factor,
    )
    actual = calculate_beta_sobolev(tau_sobolevs)
    expected = regression_data.sync_ndarray(actual)
    npt.assert_allclose(actual, expected)


def test_beta_sobolev_has_correct_thin_and_thick_limits():
    tau = pd.DataFrame(
        [[1.0e-8], [1.0e-5], [1.0e4], [1.0e6]],
        index=pd.Index(range(4), name="line"),
        columns=[0],
    )
    beta = calculate_beta_sobolev(tau)

    # Lucy's escape probability is 1 - tau/2 in the optically thin limit and
    # 1/tau in the optically thick limit, with the stable implementation's
    # transition points at 1e-4 and 1e3.
    npt.assert_allclose(beta.iloc[0, 0], 1.0 - 0.5e-8, rtol=1e-12)
    npt.assert_allclose(beta.iloc[1, 0], 1.0 - 0.5e-5, rtol=1e-12)
    npt.assert_allclose(beta.iloc[2, 0], 1.0e-4, rtol=1e-12)
    npt.assert_allclose(beta.iloc[3, 0], 1.0e-6, rtol=1e-12)
    assert np.all((beta.to_numpy() > 0.0) & (beta.to_numpy() <= 1.0))


def test_sobolev_optical_depth_scales_with_time_and_lower_population(
    nb_simulation_verysimple,
):
    legacy_plasma = nb_simulation_verysimple.plasma
    base_tau = calculate_sobolev_line_opacity(
        legacy_plasma.lines,
        legacy_plasma.level_number_density,
        legacy_plasma.time_explosion,
        legacy_plasma.stimulated_emission_factor,
    )
    doubled_tau = calculate_sobolev_line_opacity(
        legacy_plasma.lines,
        2 * legacy_plasma.level_number_density,
        2 * legacy_plasma.time_explosion,
        legacy_plasma.stimulated_emission_factor,
    )

    npt.assert_allclose(doubled_tau.to_numpy(), 4 * base_tau.to_numpy())
    assert np.all(np.isfinite(base_tau.to_numpy()))
