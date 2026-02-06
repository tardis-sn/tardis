import pytest
from tardis.visualization.tools.sdec_plot import SDECPlotter
from matplotlib.testing.compare import compare_images

# check if lineid_plot is installed
lineid_plot = pytest.importorskip("lineid_plot", reason="lineid_plot is not installed")
from tardis.visualization.tools.lineid_plotter import lineid_plotter

@pytest.fixture(scope="module")
def plotter(simulation_simple_tracked):
    """
    Create a SDECPlotter object.

    Parameters
    ----------
    simulation_simple_tracked : tardis.simulation.base.Simulation
        Simulation object.

    Returns
    -------
    tardis.visualization.tools.sdec_plot.SDECPlotter
    """
    return SDECPlotter.from_simulation(simulation_simple_tracked)


@pytest.mark.parametrize(
    "wavelengths, labels, style",
    [
        ([3951, 6355, 8567], ["Ca II", "Si II", "Ca III"], "top"),
        ([3951, 6355, 8567], ["Ca II", "Si II", "Ca III"], "inside"),
        ([3951, 6355, 8567], ["Ca II", "Si II", "Ca III"], "along"),
        pytest.param(
            [3951, 6355, 8567],  # invalid style
            ["Ca II", "Si II", "Ca III"],
            "?",
            marks=pytest.mark.xfail,
        ),
        pytest.param(  # different length arrays
            [3951, 6355],
            ["Ca II", "Si II", "Ca III"],
            "top",
            marks=pytest.mark.xfail,
        ),
    ],
)
def test_lineid_plotter(
    regression_data, plotter, tmp_path, wavelengths, labels, style
):
    ax = plotter.generate_plot_mpl()
    spectrum_wavelengths = plotter.plot_wavelength
    spectrum_data = plotter.modeled_spectrum_luminosity

    ax = lineid_plotter(
        ax,
        wavelengths,
        labels,
        spectrum_wavelengths,
        spectrum_data,
        style=style,
    )
    fig = ax.figure

    regression_data.fpath.parent.mkdir(parents=True, exist_ok=True)
    fig.figure.savefig(tmp_path / f"{regression_data.fname_prefix}.png")

    if regression_data.enable_generate_reference:
        fig.figure.savefig(
            regression_data.absolute_regression_data_dir
            / f"{regression_data.fname_prefix}.png"
        )
        pytest.skip("Skipping test to generate reference data")
    else:
        expected = str(
            regression_data.absolute_regression_data_dir
            / f"{regression_data.fname_prefix}.png"
        )
        actual = str(tmp_path / f"{regression_data.fname_prefix}.png")
        compare_images(expected, actual, tol=0.001)
