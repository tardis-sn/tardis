import logging
import os

import numpy as np
from astropy import units as u

from tardis import constants as const
from tardis.io.model.parse_geometry_configuration import (
    parse_structure_from_config,
)
from tardis.model.radiation_field_state import (
    DiluteBlackBodyRadiationFieldState,
)

logger = logging.getLogger(__name__)


def parse_radiation_field_state_from_config(
    config, geometry, dilution_factor=None, packet_source=None
):
    """
    Parses the radiation field state based on the provided configuration, radiative temperature, geometry, dilution factor, and packet source.

    Parameters
    ----------
    config : Config
        The configuration object.
    geometry : Geometry
        The geometry object.
    dilution_factor : {None, ndarray}, optional
        The dilution factor. If None, it is calculated based on the geometry.
    packet_source : {None, PacketSource}, optional
        The packet source object.

    Returns
    -------
    DiluteThermalRadiationFieldState
        The parsed radiation field state.

    Raises
    ------
    AssertionError
        If the length of t_radiative or dilution_factor is not compatible with the geometry.
    """
    (
        density_time,
        velocity,
        density,
        electron_densities,
        temperature,
    ) = parse_structure_from_config(config)

    if temperature is None:
        if config.plasma.initial_t_rad > 0 * u.K:
            temperature = (
                np.ones(geometry.no_of_shells) * config.plasma.initial_t_rad
            )
        else:
            temperature = calculate_t_radiative_from_t_inner(
                geometry, packet_source
            )

    assert len(temperature) == geometry.no_of_shells

    if dilution_factor is None:
        dilution_factor = calculate_geometric_dilution_factor(geometry)

    assert len(dilution_factor) == geometry.no_of_shells

    return DiluteBlackBodyRadiationFieldState(
        temperature, dilution_factor, geometry
    )


def parse_radiation_field_state_from_csvy(
    config, csvy_model_config, csvy_model_data, geometry, packet_source
):
    """Parses the radiation field state for CSVY model inputs.

    Parameters
    ----------
    config : Config
        The configuration object.
    csvy_model_config : Config
        CSVY model configuration.
    csvy_model_data : Config
        CSVY model data
    geometry : Geometry
        The geometry object.
    packet_source : {None, PacketSource}, optional
        The packet source object.

    Returns
    -------
    DiluteThermalRadiationFieldState
        The parsed radiation field state.
    """
    t_radiative = None
    dilution_factor = None

    if hasattr(csvy_model_data, "columns") and (
        "t_rad" in csvy_model_data.columns
    ):
        t_rad_field_index = [
            field["name"] for field in csvy_model_config.datatype.fields
        ].index("t_rad")
        t_rad_unit = u.Unit(
            csvy_model_config.datatype.fields[t_rad_field_index]["unit"]
        )
        t_radiative = csvy_model_data["t_rad"].iloc[1:].values * t_rad_unit

    elif config.plasma.initial_t_rad > 0 * u.K:
        t_radiative = (
            np.ones(geometry.no_of_shells) * config.plasma.initial_t_rad
        )
    else:
        t_radiative = calculate_t_radiative_from_t_inner(
            geometry, packet_source
        )

    if np.any(t_radiative < 1000 * u.K):
        logging.critical(
            "Radiative temperature is too low in some of the shells, temperatures below 1000K "
            f"(e.g., T_rad = {t_radiative[np.argmin(t_radiative)]} in shell {np.argmin(t_radiative)} in your model) "
            "are not accurately handled by TARDIS.",
        )

    if hasattr(csvy_model_data, "columns") and (
        "dilution_factor" in csvy_model_data.columns
    ):
        dilution_factor = csvy_model_data["dilution_factor"].iloc[1:].values
    else:
        dilution_factor = calculate_geometric_dilution_factor(geometry)

    return DiluteBlackBodyRadiationFieldState(
        t_radiative, dilution_factor, geometry
    )


def calculate_t_radiative_from_t_inner(geometry, packet_source):
    """
    Calculates the radiative temperature based on the inner temperature and the geometry of the system.

    Parameters
    ----------
    geometry : Geometry
        The geometry object.
    packet_source : PacketSource
        The packet source object.

    Returns
    -------
    Quantity
        The calculated radiative temperature.
    """
    lambda_wien_inner = const.b_wien / packet_source.temperature
    t_radiative = const.b_wien / (
        lambda_wien_inner
        * (1 + (geometry.v_middle - geometry.v_inner_boundary) / const.c)
    )
    return t_radiative


def calculate_geometric_dilution_factor(geometry):
    """Calculates the geometric dilution factor w for models without a defined w.

    Parameters
    ----------
    geometry : Geometry
        The geometry object

    Returns
    -------
    np.array
        The dilution factors for the inout geometry
    """
    return 0.5 * (
        1
        - np.sqrt(
            1
            - (
                geometry.r_inner[geometry.v_inner_boundary_index] ** 2
                / geometry.r_middle**2
            )
            .to(1)
            .value
        )
    )
