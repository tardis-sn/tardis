"""Tests for SDEC Plots."""
from tardis.base import run_tardis
import pytest
import pandas as pd
import numpy as np
import os
from copy import deepcopy
from tardis.visualization.tools.sdec_plot import SDECData, SDECPlotter
import astropy.units as u
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
import tables
import re

pytest_plugins = "regressions"

def make_valid_name(testid):
    """
    Sanitize pytest IDs to make them valid HDF group names.

    Parameters
    ----------
    testid : str
        ID to sanitize.

    Returns
    -------
    testid : str
        Sanitized ID.
    """
    testid = testid.replace("-", "_")
    testid = "_" + testid
    return testid


@pytest.fixture(scope="module")
def simulation_simple(config_verysimple, atomic_dataset):
    """
    Instantiate SDEC plotter using a simple simulation model.

    Parameters
    ----------
    config_verysimple : tardis.io.config_reader.Configuration
        Configuration object for a very simple simulation.
    atomic_dataset : str or tardis.atomic.AtomData
        Atomic data.

    Returns
    -------
    sim: tardis.simulation.base.Simulation
        Simulation object.
    """
    # Setup simulation configuration using config_verysimple and
    # override properties in such a way to make the simulation run faster
    config_verysimple.montecarlo.iterations = 3
    config_verysimple.montecarlo.no_of_packets = 4000
    config_verysimple.montecarlo.last_no_of_packets = -1
    config_verysimple.spectrum.virtual.virtual_packet_logging = True
    config_verysimple.montecarlo.no_of_virtual_packets = 1
    config_verysimple.spectrum.num = 2000
    atomic_data = deepcopy(atomic_dataset)
    sim = run_tardis(
        config_verysimple,
        atom_data=atomic_data,
        show_convergence_plots=False,
    )
    return sim


@pytest.fixture(scope="module")
def sdec_ref_data_path(tardis_ref_path):
    """
    Return the path to the reference data for the SDEC plots.

    Parameters
    ----------
    tardis_ref_path : str
        Path to the reference data directory.

    Returns
    -------
    str
        Path to SDEC reference data.
    """
    return os.path.abspath(os.path.join(tardis_ref_path, "sdec_ref.h5"))


class TestSDECPlotter:
    """Test the SDECPlotter class."""

    @pytest.fixture(scope="class")
    def plotter(self, simulation_simple):
        """
        Create a SDECPlotter object.

        Parameters
        ----------
        simulation_simple : tardis.simulation.base.Simulation
            Simulation object.

        Returns
        -------
        tardis.visualization.tools.sdec_plot.SDECPlotter
        """
        return SDECPlotter.from_simulation(simulation_simple)

    @pytest.fixture(scope="class")
    def observed_spectrum(self):
        """
        Return the observed spectrum.

        Returns
        -------
        Tuple of two astropy.units.quantity.Quantity values.
        """
        test_data = np.loadtxt(
            "tardis/visualization/tools/tests/data/observed_spectrum_test_data.dat"
        )
        observed_spectrum_wavelength, observed_spectrum_flux = test_data.T
        observed_spectrum_wavelength = observed_spectrum_wavelength * u.AA
        observed_spectrum_flux = (
            observed_spectrum_flux * u.erg / (u.s * u.cm**2 * u.AA)
        )
        return observed_spectrum_wavelength, observed_spectrum_flux

    @pytest.fixture(params=[["Si II", "Ca II", "C", "Fe I-V"]])
    def plotter_species_list(self, request, plotter):
        plotter_species_lst_obj = deepcopy(plotter)
        plotter_species_lst_obj._parse_species_list(request.param)
        yield plotter_species_lst_obj

    def test_plotter_full_species_list(
        self, plotter_species_list, data_regression, tardis_ref_path
    ):
        data_regression.check(
            plotter_species_list._full_species_list,
            fullpath=tardis_ref_path + "/" + "_full_species_list",
        )

    def test_plotter_species_list(
        self, plotter_species_list, num_regression, tardis_ref_path
    ):
        species_lst_dict = {
            "species_lst": plotter_species_list._species_list,
            # "_keep_colour": plotter_species_list._keep_colour
        }
        num_regression.check(
            species_lst_dict,
            fullpath=tardis_ref_path + "/" +"_species_list",
        )

    # @pytest.fixture(params=[["Si II", "Ca II", "C", "Fe I-V"]])
    # def plotter_species_list(self, request, plotter):
    #     plotter_calculate_plotting_data = deepcopy(plotter)
    #     plotter_calculate_plotting_data._calculate_plotting_data(request.param)
    #     yield plotter_calculate_plotting_data
    
