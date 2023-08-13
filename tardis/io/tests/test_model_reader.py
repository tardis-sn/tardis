import os

from astropy import units as u
import numpy as np
import pytest
import h5py

import tardis
from tardis.io.config_reader import Configuration
from tardis.io.model_reader import (
    read_artis_density,
    read_simple_ascii_abundances,
    read_csv_composition,
    read_uniform_abundances,
    read_cmfgen_density,
    read_cmfgen_composition,
    model_to_dict,
    transport_to_dict,
    store_model_to_hdf,
    store_transport_to_hdf,
)

data_path = os.path.join(tardis.__path__[0], "io", "tests", "data")


@pytest.fixture
def artis_density_fname():
    return os.path.join(data_path, "artis_model.dat")


@pytest.fixture
def artis_abundances_fname():
    return os.path.join(data_path, "artis_abundances.dat")


@pytest.fixture
def cmfgen_fname():
    return os.path.join(data_path, "cmfgen_model.csv")


@pytest.fixture
def csv_composition_fname():
    return os.path.join(data_path, "csv_composition.csv")


@pytest.fixture
def isotope_uniform_abundance():
    config_path = os.path.join(
        data_path, "tardis_configv1_isotope_uniabund.yml"
    )
    config = Configuration.from_yaml(config_path)
    return config.model.abundances


def test_simple_read_artis_density(artis_density_fname):
    time_of_model, velocity, mean_density = read_artis_density(
        artis_density_fname
    )

    assert np.isclose(0.00114661 * u.day, time_of_model, atol=1e-7 * u.day)
    assert np.isclose(
        mean_density[23],
        0.2250048 * u.g / u.cm**3,
        atol=1.0e-6 * u.g / u.cm**3,
    )
    assert len(mean_density) == 69
    assert len(velocity) == len(mean_density) + 1


# Artis files are currently read with read ascii files function
def test_read_simple_ascii_abundances(artis_abundances_fname):
    index, abundances = read_simple_ascii_abundances(artis_abundances_fname)
    assert len(abundances.columns) == 69
    assert np.isclose(abundances[23].loc[2], 2.672351e-08, atol=1.0e-12)


def test_read_simple_isotope_abundances(csv_composition_fname):
    index, abundances, isotope_abundance = read_csv_composition(
        csv_composition_fname
    )
    assert np.isclose(abundances.loc[6, 8], 0.5, atol=1.0e-12)
    assert np.isclose(abundances.loc[12, 5], 0.8, atol=1.0e-12)
    assert np.isclose(abundances.loc[14, 1], 0.1, atol=1.0e-12)
    assert np.isclose(isotope_abundance.loc[(28, 56), 0], 0.4, atol=1.0e-12)
    assert np.isclose(isotope_abundance.loc[(28, 58), 2], 0.7, atol=1.0e-12)
    assert abundances.shape == (4, 10)
    assert isotope_abundance.shape == (2, 10)


def test_read_cmfgen_isotope_abundances(cmfgen_fname):
    index, abundances, isotope_abundance = read_cmfgen_composition(cmfgen_fname)
    assert np.isclose(abundances.loc[6, 8], 0.5, atol=1.0e-12)
    assert np.isclose(abundances.loc[12, 5], 0.8, atol=1.0e-12)
    assert np.isclose(abundances.loc[14, 1], 0.3, atol=1.0e-12)
    assert np.isclose(isotope_abundance.loc[(28, 56), 0], 0.5, atol=1.0e-12)
    assert np.isclose(isotope_abundance.loc[(28, 58), 1], 0.7, atol=1.0e-12)
    assert abundances.shape == (4, 9)
    assert isotope_abundance.shape == (2, 9)


def test_read_uniform_abundances(isotope_uniform_abundance):
    abundances, isotope_abundance = read_uniform_abundances(
        isotope_uniform_abundance, 20
    )
    assert np.isclose(abundances.loc[8, 2], 0.19, atol=1.0e-12)
    assert np.isclose(abundances.loc[20, 5], 0.03, atol=1.0e-12)
    assert np.isclose(isotope_abundance.loc[(28, 56), 15], 0.05, atol=1.0e-12)
    assert np.isclose(isotope_abundance.loc[(28, 58), 2], 0.05, atol=1.0e-12)


def test_simple_read_cmfgen_density(cmfgen_fname):
    (
        time_of_model,
        velocity,
        mean_density,
        electron_densities,
        temperature,
    ) = read_cmfgen_density(cmfgen_fname)

    assert np.isclose(0.976 * u.day, time_of_model, atol=1e-7 * u.day)
    assert np.isclose(
        mean_density[4],
        4.2539537e-09 * u.g / u.cm**3,
        atol=1.0e-6 * u.g / u.cm**3,
    )
    assert np.isclose(
        electron_densities[5], 2.6e14 * u.cm**-3, atol=1.0e-6 * u.cm**-3
    )
    assert len(mean_density) == 9
    assert len(velocity) == len(mean_density) + 1


def test_model_to_dict(simulation_verysimple):
    model = simulation_verysimple.model

    model_dict, isotope_abundance = model_to_dict(model)

    # Check model dictionary
    assert np.array_equal(model_dict["velocity_cgs"][0], model.velocity.value)
    assert model_dict["velocity_cgs"][1] == model.velocity.unit.to_string()
    assert np.array_equal(model_dict["abundance"], model.abundance)
    assert np.array_equal(
        model_dict["time_explosion_cgs"][0], model.time_explosion.value
    )
    assert (
        model_dict["time_explosion_cgs"][1]
        == model.time_explosion.unit.to_string()
    )
    assert np.array_equal(model_dict["t_inner_cgs"][0], model.t_inner.value)
    assert model_dict["t_inner_cgs"][1] == model.t_inner.unit.to_string()
    assert np.array_equal(
        model_dict["t_radiative_cgs"][0], model.t_radiative.value
    )
    assert (
        model_dict["t_radiative_cgs"][1] == model.t_radiative.unit.to_string()
    )
    assert np.array_equal(model_dict["dilution_factor"], model.dilution_factor)
    assert np.array_equal(
        model_dict["v_boundary_inner_cgs"][0], model.v_boundary_inner.value
    )
    assert (
        model_dict["v_boundary_inner_cgs"][1]
        == model.v_boundary_inner.unit.to_string()
    )
    assert np.array_equal(
        model_dict["v_boundary_outer_cgs"][0], model.v_boundary_outer.value
    )
    assert (
        model_dict["v_boundary_outer_cgs"][1]
        == model.v_boundary_outer.unit.to_string()
    )
    assert np.array_equal(model_dict["w"], model.w)
    assert np.array_equal(model_dict["t_rad_cgs"][0], model.t_rad.value)
    assert model_dict["t_rad_cgs"][1] == model.t_rad.unit.to_string()
    assert np.array_equal(model_dict["r_inner_cgs"][0], model.r_inner.value)
    assert model_dict["r_inner_cgs"][1] == model.r_inner.unit.to_string()
    assert np.array_equal(model_dict["density_cgs"][0], model.density.value)
    assert model_dict["density_cgs"][1] == model.density.unit.to_string()


def test_store_model_to_hdf(simulation_verysimple, tmp_path):
    model = simulation_verysimple.model

    fname = tmp_path / "model.h5"

    # Store model object
    store_model_to_hdf(model, fname)

    # Check file contents
    with h5py.File(fname) as f:
        assert np.array_equal(f["model/velocity_cgs"], model.velocity.value)
        assert np.array_equal(f["model/abundance"], model.abundance)
        assert np.array_equal(
            f["model/time_explosion_cgs"], model.time_explosion.value
        )
        assert np.array_equal(f["model/t_inner_cgs"], model.t_inner.value)
        assert np.array_equal(
            f["model/t_radiative_cgs"], model.t_radiative.value
        )
        assert np.array_equal(f["model/dilution_factor"], model.dilution_factor)
        assert np.array_equal(
            f["model/v_boundary_inner_cgs"], model.v_boundary_inner.value
        )
        assert np.array_equal(
            f["model/v_boundary_outer_cgs"], model.v_boundary_outer.value
        )
        assert np.array_equal(f["model/w"], model.w)
        assert np.array_equal(f["model/t_rad_cgs"], model.t_rad.value)
        assert np.array_equal(f["model/r_inner_cgs"], model.r_inner.value)
        assert np.array_equal(f["model/density_cgs"], model.density.value)


def test_transport_to_dict(simulation_verysimple):
    transport = simulation_verysimple.transport
    transport_data = transport.__dict__

    (
        transport_dict,
        integrator_settings,
        v_packet_settings,
        virtual_spectrum_spawn_range,
    ) = transport_to_dict(transport)

    # Check transport dictionary
    for key, value in transport_dict.items():
        if isinstance(value, np.ndarray):
            if key + "_cgs" in transport_data.keys():
                assert np.array_equal(value, transport_data[key + "_cgs"])
            else:
                assert np.array_equal(value, transport_data[key])
        elif isinstance(value, list):
            assert np.array_equal(value[0], transport_data[key[:-4]].value)
            assert value[1] == transport_data[key[:-4]].unit.to_string()
        elif key == "packet_source_base_seed":  # Check packet source base seed
            assert value == transport_data["packet_source"].base_seed
        else:
            assert value == transport_data[key]

    # Check integrator settings
    for key, value in integrator_settings.items():
        assert value == transport.integrator_settings[key]

    # Check v_packet settings
    for key, value in v_packet_settings.items():
        assert value == transport.v_packet_settings[key]

    # Check virtual spectrum spawn range
    for key, value in virtual_spectrum_spawn_range.items():
        assert value.value == transport.virtual_spectrum_spawn_range[key].value


def test_store_transport_to_hdf(simulation_verysimple, tmp_path):
    transport = simulation_verysimple.transport
    transport_data = transport.__dict__
    fname = tmp_path / "transport.h5"
    # s
    # Store transport object
    store_transport_to_hdf(transport, fname)
    #
    # Check file contents
    with h5py.File(fname) as f:
        assert np.array_equal(
            f["transport/Edotlu_estimator"], transport_data["Edotlu_estimator"]
        )
        assert np.array_equal(
            f["transport/bf_heating_estimator"],
            transport_data["bf_heating_estimator"],
        )
        assert (
            f["transport/enable_full_relativity"][()]
            == transport_data["enable_full_relativity"]
        )
        assert (
            f["transport/enable_reflective_inner_boundary"][()]
            == transport_data["enable_reflective_inner_boundary"]
        )
        assert (
            f["transport/inner_boundary_albedo"][()]
            == transport_data["inner_boundary_albedo"]
        )
        assert np.array_equal(
            f["transport/input_energy"], transport_data["input_energy"]
        )
        assert np.array_equal(
            f["transport/input_mu"], transport_data["input_mu"]
        )
        assert np.array_equal(
            f["transport/input_nu"], transport_data["input_nu"]
        )
        assert np.array_equal(
            f["transport/input_r_cgs"], transport_data["input_r"].value
        )
        assert np.array_equal(
            f["transport/j_blue_estimator"], transport_data["j_blue_estimator"]
        )
        assert np.array_equal(
            f["transport/j_estimator"], transport_data["j_estimator"]
        )
        assert np.array_equal(
            f["transport/last_interaction_in_nu"],
            transport_data["last_interaction_in_nu"],
        )
        assert np.array_equal(
            f["transport/last_interaction_type"],
            transport_data["last_interaction_type"],
        )
        assert np.array_equal(
            f["transport/last_line_interaction_in_id"],
            transport_data["last_line_interaction_in_id"],
        )
        assert np.array_equal(
            f["transport/last_line_interaction_out_id"],
            transport_data["last_line_interaction_out_id"],
        )
        assert np.array_equal(
            f["transport/last_line_interaction_shell_id"],
            transport_data["last_line_interaction_shell_id"],
        )
        if hasattr(f["transport/line_interaction_type"][()], "decode"):
            assert (
                f["transport/line_interaction_type"][()].decode("utf-8")
                == transport_data["line_interaction_type"]
            )
        else:
            assert np.array_equal(
                f["transport/line_interaction_type"][()],
                transport_data["line_interaction_type"],
            )
        assert np.array_equal(
            f["transport/nu_bar_estimator"], transport_data["nu_bar_estimator"]
        )
        assert np.array_equal(
            f["transport/photo_ion_estimator"],
            transport_data["photo_ion_estimator"],
        )
        assert np.array_equal(
            f["transport/photo_ion_estimator_statistics"],
            transport_data["photo_ion_estimator_statistics"],
        )
        assert np.array_equal(
            f["transport/r_inner"], transport_data["r_inner_cgs"]
        )
        assert np.array_equal(
            f["transport/r_outer"], transport_data["r_outer_cgs"]
        )
        assert (
            f["transport/packet_source_base_seed"][()]
            == transport_data["packet_source"].base_seed
        )
        assert np.array_equal(
            f["transport/spectrum_frequency_cgs"],
            transport_data["spectrum_frequency"].value,
        )
        if hasattr(f["transport/spectrum_method"][()], "decode"):
            assert (
                f["transport/spectrum_method"][()].decode("utf-8")
                == transport_data["spectrum_method"]
            )
        else:
            assert np.array_equal(
                f["transport/spectrum_method"][()],
                transport_data["spectrum_method"],
            )
        assert np.array_equal(
            f["transport/stim_recomb_cooling_estimator"],
            transport_data["stim_recomb_cooling_estimator"],
        )
        assert np.array_equal(
            f["transport/stim_recomb_estimator"],
            transport_data["stim_recomb_estimator"],
        )
        assert np.array_equal(
            f["transport/time_of_simulation_cgs"],
            transport_data["time_of_simulation"].value,
        )
        assert f["transport/use_gpu"][()] == transport_data["use_gpu"]
        assert np.array_equal(
            f["transport/v_inner"], transport_data["v_inner_cgs"]
        )
        assert np.array_equal(
            f["transport/v_outer"], transport_data["v_outer_cgs"]
        )
        assert f["transport/nthreads"][()] == transport_data["nthreads"]
        assert f["transport/virt_logging"][()] == transport_data["virt_logging"]
        assert np.array_equal(
            f["transport/virt_packet_energies"],
            transport_data["virt_packet_energies"],
        )
        assert np.array_equal(
            f["transport/virt_packet_initial_mus"],
            transport_data["virt_packet_initial_mus"],
        )
        assert np.array_equal(
            f["transport/virt_packet_initial_rs"],
            transport_data["virt_packet_initial_rs"],
        )
        assert np.array_equal(
            f["transport/virt_packet_last_interaction_in_nu"],
            transport_data["virt_packet_last_interaction_in_nu"],
        )
        assert np.array_equal(
            f["transport/virt_packet_last_interaction_type"],
            transport_data["virt_packet_last_interaction_type"],
        )
        assert np.array_equal(
            f["transport/virt_packet_last_line_interaction_in_id"],
            transport_data["virt_packet_last_line_interaction_in_id"],
        )
        assert np.array_equal(
            f["transport/virt_packet_last_line_interaction_out_id"],
            transport_data["virt_packet_last_line_interaction_out_id"],
        )
        assert np.array_equal(
            f["transport/virt_packet_last_line_interaction_shell_id"],
            transport_data["virt_packet_last_line_interaction_shell_id"],
        )
        assert np.array_equal(
            f["transport/virt_packet_nus"], transport_data["virt_packet_nus"]
        )
        assert np.array_equal(
            f["transport/volume_cgs"], transport_data["volume"].value
        )
