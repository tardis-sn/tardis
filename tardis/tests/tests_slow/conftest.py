import os
import numpy as np
import pytest
from astropy import units as u


@pytest.fixture(scope="session")
def slow_tests_datadir():
    slow_tests_datadir = pytest.config.getvalue("slow-test-data")
    if slow_tests_datadir is None:
        pytest.skip('--slow-test-data was not specified')
    else:
        return os.path.expandvars(os.path.expanduser(slow_tests_datadir))


@pytest.fixture(scope="session")
def data_path():
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "w7")


@pytest.fixture(scope="session")
def baseline(request, slow_tests_datadir):
    """
    Fixture to ingest baseline data for slow test from already available
    compressed binaries (.npz). All data is collected in one dict and
    returned away.

    Assumed three compressed binaries (.npz) are placed in `slow-test-data/w7`
    directory, whose path is provided by command line argument:

    * ndarrays.npz                       | * quantities.npz
    Contents (all (.npy)):               | Contents (all (.npy)):
        * last_interaction_type          |     * t_rads
        * last_line_interaction_out_id   |     * luminosity_inner
        * last_line_interaction_in_id    |     * montecarlo_luminosity
        * j_estimators                   |     * montecarlo_virtual_luminousity
        * j_blue_estimators              |     * time_of_simulation
        * last_line_interaction_shell_id |     * montecarlo_nu
        * nubar_estimators               |     * last_line_interaction_angstrom
        * ws                             |     * j_blues_norm_factor

    * spectrum.npz
    Contents (all (.npy)):
        * luminosity_density_nu          | * wavelength
        * delta_frequency                | * luminosity_density_lambda
    """

    # TODO: make this fixture ingest data from an HDF5 file.
    ndarrays = dict(np.load(os.path.join(slow_tests_datadir, "ndarrays.npz")))
    quantities = dict(np.load(os.path.join(slow_tests_datadir, "quantities.npz")))
    spectrum = dict(np.load(os.path.join(slow_tests_datadir, "spectrum.npz")))

    # Associate CGS units to ndarrays of baseline quantities.
    ndarrays.update(
        j_blues_norm_factor=
            u.Quantity(quantities['j_blues_norm_factor'], '1 / (cm2 s)'),
        last_line_interaction_angstrom=
            u.Quantity(quantities['last_line_interaction_angstrom'], 'Angstrom'),
        luminosity_inner=
            u.Quantity(quantities['luminosity_inner'], 'erg / s'),
        montecarlo_luminosity=
            u.Quantity(quantities['montecarlo_luminosity'], 'erg / s'),
        montecarlo_virtual_luminosity=
            u.Quantity(quantities['montecarlo_virtual_luminosity'], 'erg / s'),
        montecarlo_nu=
            u.Quantity(quantities['montecarlo_nu'], 'Hz'),
        t_rads=
            u.Quantity(quantities['t_rads'], 'K'),
        luminosity_density_nu=
            u.Quantity(spectrum['luminosity_density_nu'], 'erg'),
        delta_frequency=
            u.Quantity(spectrum['delta_frequency'], 'Hz'),
        wavelength=
            u.Quantity(spectrum['wavelength'], 'Angstrom'),
        luminosity_density_lambda=
            u.Quantity(spectrum['luminosity_density_lambda'], 'erg / (Angstrom s)'))
    return ndarrays
