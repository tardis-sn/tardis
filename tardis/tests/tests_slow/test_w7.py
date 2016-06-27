import os
import pytest
import matplotlib.pyplot as plt
from numpy.testing import assert_allclose
from astropy.tests.helper import assert_quantity_allclose

from tardis.io.atomic import AtomData
from tardis.io.util import yaml_load_config_file
from tardis.simulation.base import Simulation
from tardis.model import Radial1DModel
from tardis.io.config_reader import Configuration


@pytest.mark.skipif(not pytest.config.getvalue("integration-tests"),
                    reason="integration tests are not included in this run")
class TestW7(object):
    """
    Slow integration test for Stratified W7 setup.
    """

    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(self, request, reference, data_path, atomic_data_fname,
              reference_datadir):
        """
        This method does initial setup of creating configuration and performing
        a single run of integration test.
        """
        self.config_file = os.path.join(data_path, "config.yml")

        # Load atom data file separately, pass it for forming tardis config.
        self.atom_data = AtomData.from_hdf5(atomic_data_fname)

        # Check whether the atom data file in current run and the atom data
        # file used in obtaining the reference data are same.
        # TODO: hard coded UUID for kurucz atom data file, generalize it later.
        kurucz_data_file_uuid1 = "5ca3035ca8b311e3bb684437e69d75d7"
        assert self.atom_data.uuid1 == kurucz_data_file_uuid1

        # Create a Configuration through yaml file and atom data.
        tardis_config = Configuration.from_yaml(
            self.config_file, atom_data=self.atom_data)

        # We now do a run with prepared config and get radial1d model.
        self.result = Radial1DModel(tardis_config)
        simulation = Simulation(tardis_config)
        simulation.legacy_run_simulation(self.result)

        # Get the reference data through the fixture.
        self.reference = reference

        # Form a base directory to save plots for `W7` setup.
        # TODO: Rough prototyping, parametrize this as more setups are added.
        self.name = "w7"
        self.base_plot_dir = os.path.join(request.config.option.tempdir, self.name)
        os.makedirs(self.base_plot_dir)

    def test_j_estimators(self):
        assert_allclose(
                self.reference['j_estimators'],
                self.result.j_estimators)

    def test_j_blue_estimators(self):
        assert_allclose(
                self.reference['j_blue_estimators'],
                self.result.j_blue_estimators)

        assert_quantity_allclose(
                self.reference['j_blues_norm_factor'],
                self.result.j_blues_norm_factor)

    def test_last_line_interactions(self):
        assert_allclose(
                self.reference['last_line_interaction_in_id'],
                self.result.last_line_interaction_in_id)

        assert_allclose(
                self.reference['last_line_interaction_out_id'],
                self.result.last_line_interaction_out_id)

        assert_allclose(
                self.reference['last_line_interaction_shell_id'],
                self.result.last_line_interaction_shell_id)

        assert_quantity_allclose(
                self.reference['last_line_interaction_angstrom'],
                self.result.last_line_interaction_angstrom)

    def test_nubar_estimators(self):
        assert_allclose(
                self.reference['nubar_estimators'],
                self.result.nubar_estimators)

    def test_ws(self):
        assert_allclose(
                self.reference['ws'],
                self.result.ws)

    def test_luminosity_inner(self):
        assert_quantity_allclose(
                self.reference['luminosity_inner'],
                self.result.luminosity_inner)

    def test_spectrum(self, plot_object):
        plot_object.add(self.plot_spectrum(), "spectrum")

        assert_quantity_allclose(
            self.reference['luminosity_density_nu'],
            self.result.runner.spectrum.luminosity_density_nu)

        assert_quantity_allclose(
            self.reference['delta_frequency'],
            self.result.runner.spectrum.delta_frequency)

        assert_quantity_allclose(
            self.reference['wavelength'],
            self.result.runner.spectrum.wavelength)

        assert_quantity_allclose(
            self.reference['luminosity_density_lambda'],
            self.result.runner.spectrum.luminosity_density_lambda)

    def plot_spectrum(self):
        plt.suptitle("Deviation in spectrum_quantities", fontweight="bold")
        figure = plt.figure()

        # `ldl_` prefixed variables associated with `luminosity_density_lambda`.
        # Axes of subplot are extracted, if we wish to make multiple plots
        # for different spectrum quantities all in one figure.
        ldl_ax = figure.add_subplot(111)
        ldl_ax.set_title("Deviation in luminosity_density_lambda")
        ldl_ax.set_xlabel("Wavelength")
        ldl_ax.set_ylabel("Relative error (1 - result / reference)")
        deviation = 1 - (
            self.result.runner.spectrum.luminosity_density_lambda.value /
            self.reference['luminosity_density_lambda'].value)

        ldl_ax.plot(self.reference['wavelength'], deviation,
                    color="blue", marker=".")

        return figure

    def test_montecarlo_properties(self):
        assert_quantity_allclose(
                self.reference['montecarlo_luminosity'],
                self.result.montecarlo_luminosity)

        assert_quantity_allclose(
                self.reference['montecarlo_virtual_luminosity'],
                self.result.runner.montecarlo_virtual_luminosity)

        assert_quantity_allclose(
                self.reference['montecarlo_nu'],
                self.result.montecarlo_nu)

    def test_shell_temperature(self, plot_object):
        plot_object.add(self.plot_t_rads(), "t_rads")

        assert_quantity_allclose(
            self.reference['t_rads'],
            self.result.t_rads)

    def plot_t_rads(self):
        plt.suptitle("Shell temperature for packets", fontweight="bold")
        figure = plt.figure()

        ax = figure.add_subplot(111)
        ax.set_xlabel("Shell id")
        ax.set_ylabel("t_rads")

        result_line = ax.plot(self.result.t_rads, color="blue",
                              marker=".", label="Result")
        reference_line = ax.plot(self.reference['t_rads'], color="green",
                                 marker=".", label="Reference")
        ax.axis([0, 28, 5000, 10000])

        error_ax = ax.twinx()
        error_line = error_ax.plot((1 - self.result.t_rads / self.reference['t_rads']),
                                   color="red", marker=".", label="Rel. Error")
        error_ax.set_ylabel("Relative error (1 - result / reference)")

        lines = result_line + reference_line + error_line
        labels = [l.get_label() for l in lines]

        ax.legend(lines, labels, loc="lower left")
        return figure
