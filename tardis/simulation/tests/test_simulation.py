import numpy.testing as npt
import pandas.util.testing as pdt
import h5py
import pytest
import os
from tardis.io.config_reader import Configuration
from tardis.model import Radial1DModel
from tardis.plasma.standard_plasmas import assemble_plasma
from tardis.simulation import Simulation
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from tardis.atomic import AtomData

@pytest.fixture
def tardis_config(tardis_config_verysimple):
    return Configuration.from_config_dict(tardis_config_verysimple)


@pytest.fixture()
def raw_model(tardis_config):
    return Radial1DModel.from_config(tardis_config)


@pytest.fixture()
def raw_plasma(tardis_config, raw_model, kurucz_atomic_data):
    return assemble_plasma(tardis_config, raw_model, kurucz_atomic_data)


@pytest.fixture()
def simulation_one_loop(raw_model, raw_plasma, tardis_config):
    sim = Simulation.from_config(tardis_config, model=raw_model,
                                 plasma=raw_plasma)
    sim.iterate(40000)

    return sim


@pytest.fixture()
def simulation_compare_data_fname():
    return 'tardis/simulation/tests/data/test_data.h5'


@pytest.fixture()
def simulation_compare_data(simulation_compare_data_fname):
    return h5py.File(simulation_compare_data_fname, mode='r')


def test_plasma_estimates(simulation_one_loop, simulation_compare_data):
    t_rad, w = simulation_one_loop.runner.calculate_radiationfield_properties()

    npt.assert_allclose(simulation_one_loop.runner.nu_bar_estimator,
                        simulation_compare_data['test1/nubar_estimators'],
                        atol=0.0)
    npt.assert_allclose(simulation_one_loop.runner.j_estimator,
                        simulation_compare_data['test1/j_estimators'],
                        atol=0.0)

    assert_quantity_allclose(
            t_rad, simulation_compare_data['test1/t_rad'] * u.Unit('K'), atol=0.0 * u.Unit('K'))
    npt.assert_allclose(w, simulation_compare_data['test1/w'], atol=0.0)


def test_packet_output(simulation_one_loop, simulation_compare_data):
    assert_quantity_allclose(
            simulation_one_loop.runner.output_nu,
            simulation_compare_data['test1/output_nu'] * u.Unit('Hz'),
            atol=0.0 * u.Unit('Hz'))

    assert_quantity_allclose(simulation_one_loop.runner.output_energy,
                        simulation_compare_data['test1/output_energy'] * u.Unit('erg'),
                        atol=0.0 * u.Unit('erg'))


def data_path(filename):
    return os.path.join('tardis/io/tests/data/', filename)


@pytest.fixture(scope="module")
def tardis_config_yml():
    filename = 'tardis_configv1_verysimple.yml'
    config = Configuration.from_yaml(data_path(filename))
    return config


@pytest.fixture(scope="module")
def atomic_data_hdf(atomic_data_fname):
    atomic_data = AtomData.from_hdf5(atomic_data_fname)

    if atomic_data.md5 != '21095dd25faa1683f4c90c911a00c3f8':
        pytest.skip('Need default Kurucz atomic dataset '
                    '(md5="21095dd25faa1683f4c90c911a00c3f8"')
    else:
        return atomic_data


@pytest.fixture(scope="module")
def simulation_config_yml(tardis_config_yml, atomic_data_hdf):
    sim = Simulation.from_config(tardis_config_yml, atom_data=atomic_data_hdf)
    sim.iterate(40000)
    return sim


@pytest.fixture(scope="module")
def hdf_file(tmpdir_factory, simulation_config_yml):
    fn = tmpdir_factory.mktemp('test_hdf').join('model_output.hdf')
    simulation_config_yml.to_hdf(str(fn))
    return fn


@pytest.fixture(scope="module")
def sim_hdf(hdf_file, atomic_data_fname):
    path = str(hdf_file)
    sim = Simulation.from_hdf(path, atomic_data_fname)
    return sim


model_quantity_attrs = ['luminosity_requested', 'time_explosion',
                        't_inner', 't_rad', 'v_inner', 'v_outer',
                        'velocity']


@pytest.mark.parametrize("attr", model_quantity_attrs)
def test_from_hdf_model_quantites(sim_hdf, simulation_config_yml, attr):
    assert_quantity_allclose(getattr(sim_hdf.model, attr), getattr(
        simulation_config_yml.model, attr))


model_nparray_attrs = ['abundance', 'dilution_factor']


@pytest.mark.parametrize("attr", model_nparray_attrs)
def test_from_hdf_model_array(sim_hdf, simulation_config_yml, attr):
    npt.assert_array_almost_equal(getattr(sim_hdf.model, attr), getattr(
        simulation_config_yml.model, attr))


plasma_properties_list = ['number_density', 'beta_rad', 'general_level_boltzmann_factor', 'level_boltzmann_factor',
                          'stimulated_emission_factor', 't_electrons', 'wavelength_cm', 'lines_lower_level_index',
                          'ionization_data', 'density', 'atomic_mass', 'level_number_density', 'lines_upper_level_index',
                          'link_t_rad_t_electron', 'nu', 'beta_sobolev', 'transition_probabilities', 'phi',
                          'electron_densities', 't_rad', 'selected_atoms', 'ion_number_density', 'partition_function',
                          'abundance', 'g_electron', 'g', 'lines',  'f_lu', 'tau_sobolevs', 'j_blues',
                          'metastability', 'w', 'time_explosion', 'excitation_energy']


@pytest.mark.parametrize("attr", plasma_properties_list)
def test_from_hdf_plasma(sim_hdf, simulation_config_yml, attr):
    if hasattr(simulation_config_yml.plasma, attr):
        npt.assert_allclose(getattr(sim_hdf.plasma, attr), getattr(
            simulation_config_yml.plasma, attr))


plasma_pd_MultiIndex = ['levels']


@pytest.mark.parametrize("attr", plasma_pd_MultiIndex)
def test_from_hdf_plasma_levels(sim_hdf, simulation_config_yml, attr):
    pdt.assert_almost_equal(getattr(sim_hdf.plasma, attr), getattr(
        simulation_config_yml.plasma, attr))


runner_properties = ['nu_bar_estimator',
                     'j_estimator',
                     'last_interaction_in_nu',
                     'last_line_interaction_in_id',
                     'last_line_interaction_out_id',
                     'last_line_interaction_shell_id',
                     'seed', 'inner_boundary_albedo',
                     'enable_reflective_inner_boundary']


@pytest.mark.parametrize("attr", runner_properties)
def test_from_hdf_runner(sim_hdf, simulation_config_yml, attr):
    if hasattr(simulation_config_yml.runner, attr):
        npt.assert_almost_equal(getattr(sim_hdf.runner, attr), getattr(
            simulation_config_yml.runner, attr))


runner_quantity_attrs = ['output_nu', 'output_energy', 'packet_luminosity',
                         'spectrum_frequency', 'sigma_thomson',
                         'montecarlo_virtual_luminosity']


@pytest.mark.parametrize("attr", runner_quantity_attrs)
def test_from_hdf_runner_quantites(sim_hdf, simulation_config_yml, attr):
    if hasattr(simulation_config_yml.runner, attr):
        assert_quantity_allclose(getattr(sim_hdf.runner, attr), getattr(
            simulation_config_yml.runner, attr))
