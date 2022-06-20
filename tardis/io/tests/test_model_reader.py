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
        if isinstance(value, np.ndarray):
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
    for key, value in integrator_settings:
        assert value == runner.integrator_settings[key]

    # Check v_packet settings
    for key, value in v_packet_settings:
        assert value == runner.v_packet_settings[key]

    # Check virtual spectrum spawn range
    for key, value in virtual_spectrum_spawn_range:
        assert value.value == runner.virtual_spectrum_spawn_range[key].value
