# reading different model files

from tardis.io.configuration.config_reader import ConfigurationNameSpace
from tardis.montecarlo.base import MontecarloTransport
from tardis.montecarlo.packet_source import (
    BlackBodySimpleSource,
    BlackBodySimpleSourceRelativistic,
)
from tardis.io.model.density import calculate_density_after_time
from astropy import units as u
import h5py

import logging

# Adding logging support
logger = logging.getLogger(__name__)


def transport_to_dict(transport):
    """
    Retrieves all the data from a transport object and returns a dictionary.

    Parameters
    ----------
    transport : tardis.montecarlo.MontecarloTransport

    Returns
    -------
    transport_dict : dict
    integrator_settings : dict
    v_packet_settings : dict
    virtual_spectrum_spawn_range : dict
    """

    transport_dict = {
        "Edotlu_estimator": transport.Edotlu_estimator,
        "bf_heating_estimator": transport.bf_heating_estimator,
        "disable_electron_scattering": transport.disable_electron_scattering,
        "enable_full_relativity": transport.enable_full_relativity,
        "enable_reflective_inner_boundary": transport.enable_reflective_inner_boundary,
        "inner_boundary_albedo": transport.inner_boundary_albedo,
        "input_energy": transport.input_energy,
        "input_mu": transport.input_mu,
        "input_nu": transport.input_nu,
        "input_r": transport.input_r,
        "j_blue_estimator": transport.j_blue_estimator,
        "j_estimator": transport.j_estimator,
        "last_interaction_in_nu": transport.last_interaction_in_nu,
        "last_interaction_type": transport.last_interaction_type,
        "last_line_interaction_in_id": transport.last_line_interaction_in_id,
        "last_line_interaction_out_id": transport.last_line_interaction_out_id,
        "last_line_interaction_shell_id": transport.last_line_interaction_shell_id,
        "line_interaction_type": transport.line_interaction_type,
        "nu_bar_estimator": transport.nu_bar_estimator,
        "photo_ion_estimator": transport.photo_ion_estimator,
        "photo_ion_estimator_statistics": transport.photo_ion_estimator_statistics,
        "r_inner": transport.r_inner_cgs,
        "r_outer": transport.r_outer_cgs,
        "packet_source_base_seed": transport.packet_source.base_seed,
        "spectrum_frequency_cgs": transport.spectrum_frequency,
        "spectrum_method": transport.spectrum_method,
        "stim_recomb_cooling_estimator": transport.stim_recomb_cooling_estimator,
        "stim_recomb_estimator": transport.stim_recomb_estimator,
        "time_of_simulation_cgs": transport.time_of_simulation,
        "use_gpu": transport.use_gpu,
        "v_inner": transport.v_inner_cgs,
        "v_outer": transport.v_outer_cgs,
        "nthreads": transport.nthreads,
        "virt_logging": transport.virt_logging,
        "virt_packet_energies": transport.virt_packet_energies,
        "virt_packet_initial_mus": transport.virt_packet_initial_mus,
        "virt_packet_initial_rs": transport.virt_packet_initial_rs,
        "virt_packet_last_interaction_in_nu": transport.virt_packet_last_interaction_in_nu,
        "virt_packet_last_interaction_type": transport.virt_packet_last_interaction_type,
        "virt_packet_last_line_interaction_in_id": transport.virt_packet_last_line_interaction_in_id,
        "virt_packet_last_line_interaction_out_id": transport.virt_packet_last_line_interaction_out_id,
        "virt_packet_last_line_interaction_shell_id": transport.virt_packet_last_line_interaction_shell_id,
        "virt_packet_nus": transport.virt_packet_nus,
        "volume_cgs": transport.volume,
    }

    for key, value in transport_dict.items():
        if key.endswith("_cgs"):
            transport_dict[key] = [value.cgs.value, value.unit.to_string()]

    integrator_settings = transport.integrator_settings
    v_packet_settings = transport.v_packet_settings
    virtual_spectrum_spawn_range = transport.virtual_spectrum_spawn_range

    return (
        transport_dict,
        integrator_settings,
        v_packet_settings,
        virtual_spectrum_spawn_range,
    )


def store_transport_to_hdf(transport, fname):
    """
    Stores data from MontecarloTransport object into a hdf file.

    Parameters
    ----------
    transport : tardis.montecarlo.MontecarloTransport
    filename : str
    """

    with h5py.File(fname, "a") as f:
        transport_group = f.require_group("transport")
        transport_group.clear()

        (
            transport_data,
            integrator_settings,
            v_packet_settings,
            virtual_spectrum_spawn_range,
        ) = transport_to_dict(transport)

        for key, value in transport_data.items():
            if key.endswith("_cgs"):
                transport_group.create_dataset(key, data=value[0])
                transport_group.create_dataset(key + "_unit", data=value[1])
            else:
                if value is not None:
                    transport_group.create_dataset(key, data=value)

        integrator_settings_group = transport_group.create_group(
            "integrator_settings"
        )
        for key, value in integrator_settings.items():
            integrator_settings_group.create_dataset(key, data=value)

        v_packet_settings_group = transport_group.create_group(
            "v_packet_settings"
        )
        for key, value in v_packet_settings.items():
            v_packet_settings_group.create_dataset(key, data=value)

        virtual_spectrum_spawn_range_group = transport_group.create_group(
            "virtual_spectrum_spawn_range"
        )
        for key, value in virtual_spectrum_spawn_range.items():
            virtual_spectrum_spawn_range_group.create_dataset(key, data=value)


def transport_from_hdf(fname):
    """
    Creates a MontecarloTransport object using data stored in a hdf file.

    Parameters
    ----------
    fname : str

    Returns
    -------
    new_transport : tardis.montecarlo.MontecarloTransport
    """

    d = {}

    # Loading data from hdf file
    with h5py.File(fname, "r") as f:
        transport_group = f["transport"]
        for key, value in transport_group.items():
            if not key.endswith("_unit"):
                if isinstance(value, h5py.Dataset):
                    if isinstance(value[()], bytes):
                        d[key] = value[()].decode("utf-8")
                    else:
                        d[key] = value[()]
                else:
                    data_inner = {}
                    for key_inner, value_inner in value.items():
                        if isinstance(value_inner[()], bytes):
                            data_inner[key] = value_inner[()].decode("utf-8")
                        else:
                            data_inner[key] = value_inner[()]
                        data_inner[key_inner] = value_inner[()]
                    d[key] = data_inner

        for key, value in transport_group.items():
            if key.endswith("_unit"):
                d[key[:-5]] = [d[key[:-5]], value[()]]

    # Converting cgs data to astropy quantities
    for key, value in d.items():
        if key.endswith("_cgs"):
            d[key] = u.Quantity(value[0], unit=u.Unit(value[1].decode("utf-8")))

    # Using packet source seed to packet source
    if not d["enable_full_relativity"]:
        d["packet_source"] = BlackBodySimpleSource(d["packet_source_base_seed"])
    else:
        d["packet_source"] = BlackBodySimpleSourceRelativistic(
            d["packet_source_base_seed"]
        )

    # Converting virtual spectrum spawn range values to astropy quantities
    vssr = d["virtual_spectrum_spawn_range"]
    d["virtual_spectrum_spawn_range"] = {
        "start": u.Quantity(vssr["start"], unit=u.Unit("Angstrom")),
        "end": u.Quantity(vssr["end"], unit=u.Unit("Angstrom")),
    }

    # Converting dictionaries to ConfigurationNameSpace
    d["integrator_settings"] = ConfigurationNameSpace(d["integrator_settings"])
    d["v_packet_settings"] = ConfigurationNameSpace(d["v_packet_settings"])
    d["virtual_spectrum_spawn_range"] = ConfigurationNameSpace(
        d["virtual_spectrum_spawn_range"]
    )

    # Creating a transport object and storing data
    new_transport = MontecarloTransport(
        spectrum_frequency=d["spectrum_frequency_cgs"],
        virtual_spectrum_spawn_range=d["virtual_spectrum_spawn_range"],
        disable_electron_scattering=d["disable_electron_scattering"],
        enable_reflective_inner_boundary=d["enable_reflective_inner_boundary"],
        enable_full_relativity=d["enable_full_relativity"],
        inner_boundary_albedo=d["inner_boundary_albedo"],
        line_interaction_type=d["line_interaction_type"],
        integrator_settings=d["integrator_settings"],
        v_packet_settings=d["v_packet_settings"],
        spectrum_method=d["spectrum_method"],
        packet_source=d["packet_source"],
        nthreads=d["nthreads"],
        virtual_packet_logging=d["virt_logging"],
        use_gpu=d["use_gpu"],
    )

    new_transport.Edotlu_estimator = d["Edotlu_estimator"]
    new_transport.bf_heating_estimator = d["bf_heating_estimator"]
    new_transport.input_energy = d["input_energy"]
    new_transport.input_mu = d["input_mu"]
    new_transport.input_nu = d["input_nu"]
    new_transport.input_r = d["input_r_cgs"]
    new_transport.j_blue_estimator = d["j_blue_estimator"]
    new_transport.j_estimator = d["j_estimator"]
    new_transport.last_interaction_in_nu = d["last_interaction_in_nu"]
    new_transport.last_interaction_type = d["last_interaction_type"]
    new_transport.last_line_interaction_in_id = d["last_line_interaction_in_id"]
    new_transport.last_line_interaction_out_id = d[
        "last_line_interaction_out_id"
    ]
    new_transport.last_line_interaction_shell_id = d[
        "last_line_interaction_shell_id"
    ]
    new_transport.nu_bar_estimator = d["nu_bar_estimator"]
    new_transport.photo_ion_estimator = d["photo_ion_estimator"]
    new_transport.photo_ion_estimator_statistics = d[
        "photo_ion_estimator_statistics"
    ]
    new_transport.r_inner_cgs = d["r_inner"]
    new_transport.r_outer_cgs = d["r_outer"]
    new_transport.stim_recomb_cooling_estimator = d[
        "stim_recomb_cooling_estimator"
    ]
    new_transport.stim_recomb_estimator = d["stim_recomb_estimator"]
    new_transport.time_of_simulation = d["time_of_simulation_cgs"]
    new_transport.v_inner_cgs = d["v_inner"]
    new_transport.v_outer_cgs = d["v_outer"]
    new_transport.virt_packet_energies = d["virt_packet_energies"]
    new_transport.virt_packet_initial_mus = d["virt_packet_initial_mus"]
    new_transport.virt_packet_initial_rs = d["virt_packet_initial_rs"]
    new_transport.virt_packet_last_interaction_in_nu = d[
        "virt_packet_last_interaction_in_nu"
    ]
    new_transport.virt_packet_last_interaction_type = d[
        "virt_packet_last_interaction_type"
    ]
    new_transport.virt_packet_last_line_interaction_in_id = d[
        "virt_packet_last_line_interaction_in_id"
    ]
    new_transport.virt_packet_last_line_interaction_out_id = d[
        "virt_packet_last_line_interaction_out_id"
    ]
    new_transport.virt_packet_last_line_interaction_shell_id = d[
        "virt_packet_last_line_interaction_shell_id"
    ]
    new_transport.virt_packet_nus = d["virt_packet_nus"]
    new_transport.volume = d["volume_cgs"]

    return new_transport


def simulation_state_to_dict(simulation_state):
    """
    Retrieves all the data from a SimulationState object and returns a dictionary.

    Parameters
    ----------
    transport : tardis.model.SimulationState

    Returns
    -------
    model_dict : dict
    isotope_abundance : dict
    """
    simulation_state_dict = {
        "velocity_cgs": simulation_state.velocity.cgs,
        "abundance": simulation_state.abundance,
        "time_explosion_cgs": simulation_state.time_explosion.cgs,
        "t_inner_cgs": simulation_state.t_inner.cgs,
        "t_radiative_cgs": simulation_state.t_radiative.cgs,
        "dilution_factor": simulation_state.dilution_factor,
        "v_boundary_inner_cgs": simulation_state.v_boundary_inner.cgs,
        "v_boundary_outer_cgs": simulation_state.v_boundary_outer.cgs,
        "dilution_factor": simulation_state.dilution_factor,
        "r_inner_cgs": simulation_state.r_inner.cgs,
        "density_cgs": simulation_state.density.cgs,
    }

    for key, value in simulation_state_dict.items():
    for key, value in simulation_state_dict.items():
        if hasattr(value, "unit"):
            simulation_state_dict[key] = [
                value.cgs.value,
                value.cgs.unit.to_string(),
            ]

    return simulation_state_dict
