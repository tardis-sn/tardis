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
    runner_to_dict,
    store_model_to_hdf,
    store_runner_to_hdf,
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

    model_dict, homologous_density, isotope_abundance = model_to_dict(model)

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

    # Check homologous density
    assert np.array_equal(
        homologous_density["optional_hdf_properties"],
        model.homologous_density.optional_hdf_properties,
    )
    assert np.array_equal(
        homologous_density["density_0"], model.homologous_density.density_0
    )
    assert np.array_equal(
        homologous_density["time_0"], model.homologous_density.time_0
    )


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

        # Check homologous density
        assert np.array_equal(
            f["model/homologous_density/optional_hdf_properties"],
            model.homologous_density.optional_hdf_properties,
        )
        assert np.array_equal(
            f["model/homologous_density/density_0"],
            model.homologous_density.density_0.value,
        )
        assert np.array_equal(
            f["model/homologous_density/time_0"],
            model.homologous_density.time_0.value,
        )


def test_runner_to_dict(simulation_verysimple):
    runner = simulation_verysimple.runner
    runner_data = runner.__dict__

    (
        runner_dict,
        integrator_settings,
        v_packet_settings,
        virtual_spectrum_spawn_range,
    ) = runner_to_dict(runner)

    # Check runner dictionary
    for key, value in runner_dict.items():
        if key == "single_packet_seed":
            if value is None:
                assert key not in runner_data.keys()
            else:
                assert value == runner_data[key]
        elif isinstance(value, np.ndarray):
            if key + "_cgs" in runner_data.keys():
                assert np.array_equal(value, runner_data[key + "_cgs"])
            else:
                assert np.array_equal(value, runner_data[key])
        elif isinstance(value, list):
            assert np.array_equal(value[0], runner_data[key[:-4]].value)
            assert value[1] == runner_data[key[:-4]].unit.to_string()
        else:
            assert value == runner_data[key]

    # Check integrator settings
    for key, value in integrator_settings.items():
        assert value == runner.integrator_settings[key]

    # Check v_packet settings
    for key, value in v_packet_settings.items():
        assert value == runner.v_packet_settings[key]

    # Check virtual spectrum spawn range
    for key, value in virtual_spectrum_spawn_range.items():
        assert value.value == runner.virtual_spectrum_spawn_range[key].value


def test_store_runner_to_hdf(simulation_verysimple, tmp_path):
    runner = simulation_verysimple.runner
    runner_data = runner.__dict__
    fname = tmp_path / "runner.h5"
    # s
    # Store runner object
    store_runner_to_hdf(runner, fname)
    #
    # Check file contents
    with h5py.File(fname) as f:
        assert np.array_equal(
            f["runner/Edotlu_estimator"], runner_data["Edotlu_estimator"]
        )
        assert np.array_equal(
            f["runner/bf_heating_estimator"],
            runner_data["bf_heating_estimator"],
        )
        assert (
            f["runner/enable_full_relativity"][()]
            == runner_data["enable_full_relativity"]
        )
        assert (
            f["runner/enable_reflective_inner_boundary"][()]
            == runner_data["enable_reflective_inner_boundary"]
        )
        assert (
            f["runner/inner_boundary_albedo"][()]
            == runner_data["inner_boundary_albedo"]
        )
        assert np.array_equal(
            f["runner/input_energy"], runner_data["input_energy"]
        )
        assert np.array_equal(f["runner/input_mu"], runner_data["input_mu"])
        assert np.array_equal(f["runner/input_nu"], runner_data["input_nu"])
        assert np.array_equal(
            f["runner/input_r_cgs"], runner_data["input_r"].value
        )
        assert np.array_equal(
            f["runner/j_blue_estimator"], runner_data["j_blue_estimator"]
        )
        assert np.array_equal(
            f["runner/j_estimator"], runner_data["j_estimator"]
        )
        assert np.array_equal(
            f["runner/last_interaction_in_nu"],
            runner_data["last_interaction_in_nu"],
        )
        assert np.array_equal(
            f["runner/last_interaction_type"],
            runner_data["last_interaction_type"],
        )
        assert np.array_equal(
            f["runner/last_line_interaction_in_id"],
            runner_data["last_line_interaction_in_id"],
        )
        assert np.array_equal(
            f["runner/last_line_interaction_out_id"],
            runner_data["last_line_interaction_out_id"],
        )
        assert np.array_equal(
            f["runner/last_line_interaction_shell_id"],
            runner_data["last_line_interaction_shell_id"],
        )
        if hasattr(f["runner/line_interaction_type"][()], "decode"):
            assert (
                f["runner/line_interaction_type"][()].decode("utf-8")
                == runner_data["line_interaction_type"]
            )
        else:
            assert np.array_equal(
                f["runner/line_interaction_type"][()],
                runner_data["line_interaction_type"],
            )
        assert np.array_equal(
            f["runner/nu_bar_estimator"], runner_data["nu_bar_estimator"]
        )
        assert np.array_equal(
            f["runner/photo_ion_estimator"], runner_data["photo_ion_estimator"]
        )
        assert np.array_equal(
            f["runner/photo_ion_estimator_statistics"],
            runner_data["photo_ion_estimator_statistics"],
        )
        assert np.array_equal(f["runner/r_inner"], runner_data["r_inner_cgs"])
        assert np.array_equal(f["runner/r_outer"], runner_data["r_outer_cgs"])
        assert f["runner/seed"][()] == runner_data["seed"]
        assert np.array_equal(
            f["runner/spectrum_frequency_cgs"],
            runner_data["spectrum_frequency"].value,
        )
        if hasattr(f["runner/spectrum_method"][()], "decode"):
            assert (
                f["runner/spectrum_method"][()].decode("utf-8")
                == runner_data["spectrum_method"]
            )
        else:
            assert np.array_equal(
                f["runner/spectrum_method"][()],
                runner_data["spectrum_method"],
            )
        assert np.array_equal(
            f["runner/stim_recomb_cooling_estimator"],
            runner_data["stim_recomb_cooling_estimator"],
        )
        assert np.array_equal(
            f["runner/stim_recomb_estimator"],
            runner_data["stim_recomb_estimator"],
        )
        assert np.array_equal(
            f["runner/time_of_simulation_cgs"],
            runner_data["time_of_simulation"].value,
        )
        assert f["runner/use_gpu"][()] == runner_data["use_gpu"]
        assert np.array_equal(f["runner/v_inner"], runner_data["v_inner_cgs"])
        assert np.array_equal(f["runner/v_outer"], runner_data["v_outer_cgs"])
        assert f["runner/virt_logging"][()] == runner_data["virt_logging"]
        assert np.array_equal(
            f["runner/virt_packet_energies"],
            runner_data["virt_packet_energies"],
        )
        assert np.array_equal(
            f["runner/virt_packet_initial_mus"],
            runner_data["virt_packet_initial_mus"],
        )
        assert np.array_equal(
            f["runner/virt_packet_initial_rs"],
            runner_data["virt_packet_initial_rs"],
        )
        assert np.array_equal(
            f["runner/virt_packet_last_interaction_in_nu"],
            runner_data["virt_packet_last_interaction_in_nu"],
        )
        assert np.array_equal(
            f["runner/virt_packet_last_interaction_type"],
            runner_data["virt_packet_last_interaction_type"],
        )
        assert np.array_equal(
            f["runner/virt_packet_last_line_interaction_in_id"],
            runner_data["virt_packet_last_line_interaction_in_id"],
        )
        assert np.array_equal(
            f["runner/virt_packet_last_line_interaction_out_id"],
            runner_data["virt_packet_last_line_interaction_out_id"],
        )
        assert np.array_equal(
            f["runner/virt_packet_nus"], runner_data["virt_packet_nus"]
        )
        assert np.array_equal(
            f["runner/volume_cgs"], runner_data["volume"].value
        )
        if "runner/single_packet_seed" in f:
            assert (
                f["runner/single_packet_seed"][()]
                == runner_data["single_packet_seed"]
            )
