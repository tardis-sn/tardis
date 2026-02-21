import logging

import numpy as np
import pandas as pd
from astropy import units as u

from tardis import constants as const
from tardis.io.configuration.config_reader import (
    Configuration,
)
from tardis.io.model.parse_geometry_configuration import (
    parse_structure_from_config,
)
from tardis.model.geometry.radial1d import HomologousRadial1DGeometry
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.radiation_field.validate_radiation_field import (
    validate_radiative_temperature,
)
from tardis.transport.montecarlo.packet_source.base import BasePacketSource

logger = logging.getLogger(__name__)


def parse_radiation_field_state_from_config(
    config: Configuration,
    geometry: HomologousRadial1DGeometry,
    dilution_factor: np.ndarray | None = None,
    packet_source: BasePacketSource | None = None,
) -> DilutePlanckianRadiationField:
    """Parses the radiation field state based on the provided configuration,
    geometry, dilution factor, and packet source.

    Parameters
    ----------
    config
        The configuration object.
    geometry
        The geometry object.
    dilution_factor
        The dilution factor. If None, it is calculated based on the geometry.
    packet_source
        The packet source object.

    Returns
    -------
    radiation_field
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

    return DilutePlanckianRadiationField(temperature, dilution_factor, geometry)


def parse_radiation_field_state_from_csvy(
    config: Configuration,
    csvy_model_config: Configuration,
    csvy_model_data: pd.DataFrame | None,
    geometry: HomologousRadial1DGeometry,
    packet_source: BasePacketSource,
) -> DilutePlanckianRadiationField:
    """Parses the radiation field state for CSVY model inputs.

    Parameters
    ----------
    config
        The configuration object.
    csvy_model_config
        CSVY model configuration.
    csvy_model_data
        CSVY model data.
    geometry
        The geometry object.
    packet_source
        The packet source object.

    Returns
    -------
    radiation_field
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

    validate_radiative_temperature(t_radiative)

    if hasattr(csvy_model_data, "columns") and (
        "dilution_factor" in csvy_model_data.columns
    ):
        dilution_factor = csvy_model_data["dilution_factor"].iloc[1:].values
    else:
        dilution_factor = calculate_geometric_dilution_factor(geometry)

    return DilutePlanckianRadiationField(t_radiative, dilution_factor, geometry)


def calculate_t_radiative_from_t_inner(
    geometry: HomologousRadial1DGeometry,
    packet_source: BasePacketSource,
) -> u.Quantity:
    """Calculates the radiative temperature based on the inner temperature
    and the geometry of the system.

    Parameters
    ----------
    geometry
        The geometry object.
    packet_source
        The packet source object.

    Returns
    -------
    t_radiative
        The calculated radiative temperature.
    """
    lambda_wien_inner = const.b_wien / packet_source.temperature
    t_radiative = const.b_wien / (
        lambda_wien_inner
        * (1 + (geometry.v_middle - geometry.v_inner_boundary) / const.c)
    )
    return t_radiative


def calculate_geometric_dilution_factor(
    geometry: HomologousRadial1DGeometry,
) -> np.ndarray:
    """Calculates the geometric dilution factor w for models without a defined w.

    Parameters
    ----------
    geometry
        The geometry object.

    Returns
    -------
    dilution_factors
        The dilution factors for the input geometry.
    """
    value = (
        1 - (geometry.r_inner_active[0] ** 2 / geometry.r_middle**2).to(1).value
    )
    value = np.clip(value, 0, None)
    return 0.5 * (1 - np.sqrt(value))
