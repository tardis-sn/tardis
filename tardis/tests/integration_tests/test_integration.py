import os

import yaml
import pytest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpy.testing import assert_allclose

from tardis.io.atom_data.base import AtomData
from tardis.simulation import Simulation
from tardis.io.config_reader import Configuration

quantity_comparison = [
    ('/simulation/runner/last_line_interaction_in_id',
     'runner.last_line_interaction_in_id'),
    ('/simulation/runner/last_line_interaction_out_id',
     'runner.last_line_interaction_out_id'),
    ('/simulation/runner/last_line_interaction_shell_id',
     'runner.last_line_interaction_shell_id'),
    ('/simulation/plasma/j_blues',
     'plasma.j_blues'),
    ('/simulation/plasma/j_blue_estimator',
     'plasma.j_blue_estimator'),
    ('/simulation/runner/packet_luminosity',
     'runner.packet_luminosity.cgs.value'),
    ('/simulation/runner/montecarlo_virtual_luminosity',
     'runner.montecarlo_virtual_luminosity.cgs.value'),
    ('/simulation/runner/output_nu',
     'runner.output_nu.cgs.value'),
    ('/simulation/plasma/ion_number_density',
     'plasma.ion_number_density'),
    ('/simulation/plasma/level_number_density',
     'plasma.level_number_density'),
    ('/simulation/plasma/electron_densities',
     'plasma.electron_densities'),
    ('/simulation/plasma/tau_sobolevs',
     'plasma.tau_sobolevs'),
    ('/simulation/plasma/transition_probabilities',
     'plasma.transition_probabilities'),
    ('/simulation/model/t_radiative',
     'model.t_radiative.cgs.value'),
    ('/simulation/model/w',
     'model.w'),
    ('/simulation/runner/j_estimator',
     'runner.j_estimator'),
    ('/simulation/runner/nu_bar_estimator',
     'runner.nu_bar_estimator'),
    ('/simulation/plasma/j_blues_norm_factor',
     'plasma.j_blues_norm_factor.cgs.value'),
    ('/simulation/plasma/luminosity_inner',
     'plasma.luminosity_inner.cgs.value'),
     ]


@pytest.fixture(params=quantity_comparison)
def model_quantities(request):
    return request.param


@pytest.mark.skipif(not pytest.config.getvalue("integration-tests"),
                    reason="integration tests are not included in this run")
@pytest.mark.integration
class TestIntegration(object):
    """Slow integration test for various setups present in subdirectories of
    ``tardis/tests/integration_tests``.
    """

    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(self, request, reference, data_path, pytestconfig):
        """
        This method does initial setup of creating configuration and performing
        a single run of integration test.
        """
        # Get capture manager
        capmanager = pytestconfig.pluginmanager.getplugin('capturemanager')

        # The last component in dirpath can be extracted as name of setup.
        self.name = data_path['setup_name']

        self.config_file = os.path.join(data_path['config_dirpath'], "config.yml")

        # A quick hack to use atom data per setup. Atom data is ingested from
        # local HDF or downloaded and cached from a url, depending on data_path
        # keys.
        atom_data_name = yaml.load(open(self.config_file))['atom_data']

        # Get the path to HDF file:
        atom_data_filepath = os.path.join(
            data_path['atom_data_path'], atom_data_name
        )

        # Load atom data file separately, pass it for forming tardis config.
        self.atom_data = AtomData.from_hdf(atom_data_filepath)

        # Check whether the atom data file in current run and the atom data
        # file used in obtaining the reference data are same.
        # TODO: hard coded UUID for kurucz atom data file, generalize it later.
        # kurucz_data_file_uuid1 = "5ca3035ca8b311e3bb684437e69d75d7"
        # assert self.atom_data.uuid1 == kurucz_data_file_uuid1

        # Create a Configuration through yaml file and atom data.
        tardis_config = Configuration.from_yaml(self.config_file)

        # Check whether current run is with less packets.
        if request.config.getoption("--less-packets"):
            less_packets = request.config.integration_tests_config['less_packets']
            tardis_config['montecarlo']['no_of_packets'] = (
                less_packets['no_of_packets']
            )
            tardis_config['montecarlo']['last_no_of_packets'] = (
                less_packets['last_no_of_packets']
            )




        # We now do a run with prepared config and get the simulation object.
        self.result = Simulation.from_config(tardis_config,
                                             atom_data=self.atom_data)

        capmanager.suspendcapture(True)
        # If current test run is just for collecting reference data, store the
        # output model to HDF file, save it at specified path. Skip all tests.
        # Else simply perform the run and move further for performing
        # assertions.
        self.result.run()
        if request.config.getoption("--generate-reference"):
            ref_data_path = os.path.join(
                data_path['reference_path'], "{0}.h5".format(self.name)
            )
            if os.path.exists(ref_data_path):
                pytest.skip(
                    'Reference data {0} does exist and tests will not '
                    'proceed generating new data'.format(ref_data_path))
            self.result.to_hdf(file_path=ref_data_path)
            pytest.skip("Reference data saved at {0}".format(
                data_path['reference_path']
            ))
        capmanager.resumecapture()

        # Get the reference data through the fixture.
        self.reference = reference

    def test_model_quantities(self, model_quantities):
        reference_quantity_name, tardis_quantity_name = model_quantities
        if reference_quantity_name not in self.reference:
            pytest.skip('{0} not calculated in this run'.format(
                reference_quantity_name))
        reference_quantity = self.reference[reference_quantity_name]
        tardis_quantity = eval('self.result.' + tardis_quantity_name)
        assert_allclose(tardis_quantity, reference_quantity)

    def plot_t_rad(self):
        plt.suptitle("Shell temperature for packets", fontweight="bold")
        figure = plt.figure()

        ax = figure.add_subplot(111)
        ax.set_xlabel("Shell id")
        ax.set_ylabel("t_rad")

        result_line = ax.plot(
            self.result.model.t_rad.cgs, color="blue", marker=".", label="Result"
        )
        reference_line = ax.plot(
            self.reference['/simulation/model/t_rad'],
            color="green", marker=".", label="Reference"
        )

        error_ax = ax.twinx()
        error_line = error_ax.plot(
            (1 - self.result.model.t_rad.cgs.value / self.reference['/simulation/model/t_rad']),
            color="red", marker=".", label="Rel. Error"
        )
        error_ax.set_ylabel("Relative error (1 - result / reference)")

        lines = result_line + reference_line + error_line
        labels = [l.get_label() for l in lines]

        ax.legend(lines, labels, loc="lower left")
        return figure


    def test_spectrum(self, plot_object):
        plot_object.add(self.plot_spectrum(), "{0}_spectrum".format(self.name))

        assert_allclose(
            self.reference['/simulation/runner/spectrum/luminosity_density_nu'],
            self.result.runner.spectrum.luminosity_density_nu.cgs.value)

        assert_allclose(
            self.reference['/simulation/runner/spectrum/wavelength'],
            self.result.runner.spectrum.wavelength.cgs.value)

        assert_allclose(
            self.reference['/simulation/runner/spectrum/luminosity_density_lambda'],
            self.result.runner.spectrum.luminosity_density_lambda.cgs.value)

    def plot_spectrum(self):

        # `ldl_` prefixed variables associated with `luminosity_density_lambda`.
        # Axes of subplot are extracted, if we wish to make multiple plots
        # for different spectrum quantities all in one figure.
        gs = plt.GridSpec(2, 1, height_ratios=[3, 1])

        spectrum_ax = plt.subplot(gs[0])

        spectrum_ax.set_ylabel("Flux [cgs]")
        deviation = 1 - (
            self.result.runner.spectrum.luminosity_density_lambda.cgs.value /
            self.reference[
                '/simulation/runner/spectrum/luminosity_density_lambda']

        )


        spectrum_ax.plot(
            self.reference['/simulation/runner/spectrum/wavelength'],
            self.reference[
                '/simulation/runner/spectrum/luminosity_density_lambda'],
            color="black"
        )

        spectrum_ax.plot(
            self.reference['/simulation/runner/spectrum/wavelength'],
            self.result.runner.spectrum.luminosity_density_lambda.cgs.value,
            color="red"
        )
        spectrum_ax.set_xticks([])
        deviation_ax = plt.subplot(gs[1])
        deviation_ax.plot(self.reference['/simulation/runner/spectrum/wavelength'],
                   deviation, color='black')
        deviation_ax.set_xlabel("Wavelength [Angstrom]")

        return plt.gcf()
