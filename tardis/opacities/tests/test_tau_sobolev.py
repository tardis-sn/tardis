import numpy.testing as npt
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
