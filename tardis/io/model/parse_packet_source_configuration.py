from astropy import units as u

from tardis.transport.montecarlo.packet_source import (
    BlackBodySimpleSource,
    BlackBodySimpleSourceRelativistic,
)


def initialize_packet_source(packet_source, config, geometry):
    """
    Initialize the packet source based on config and geometry

    Parameters
    ----------
    config : Config
        The configuration object containing the supernova and plasma settings.
    geometry : Geometry
        The geometry object containing the inner radius information.
    packet_source : BasePacketSource
        The packet source object based on the configuration and geometry.

    Returns
    -------
    packet_source : BasePacketSource
        The packet source object based on the configuration and geometry.

    Raises
    ------
    ValueError
        If both t_inner and luminosity_requested are None.
    """
    luminosity_requested = config.supernova.luminosity_requested
    if config.plasma.initial_t_inner > 0.0 * u.K:
        packet_source.radius = geometry.r_inner_active[0]
        packet_source.temperature = config.plasma.initial_t_inner

    elif (config.plasma.initial_t_inner < 0.0 * u.K) and (
        luminosity_requested is not None
    ):
        packet_source.radius = geometry.r_inner_active[0]
        packet_source.set_temperature_from_luminosity(luminosity_requested)
    else:
        raise ValueError(
            "Both t_inner and luminosity_requested cannot be None."
        )
    return packet_source


def parse_packet_source_from_config(config, geometry, legacy_mode_enabled):
    """
    Parse the packet source based on the given configuration and geometry.

    Parameters
    ----------
    config : Config
        The configuration object containing the supernova and plasma settings.
    geometry : Geometry
        The geometry object containing the inner radius information.

    Returns
    -------
    packet_source : BlackBodySimpleSource
        The packet source object based on the configuration and geometry.
    """
    if config.montecarlo.enable_full_relativity:
        packet_source = BlackBodySimpleSourceRelativistic(
            base_seed=config.montecarlo.seed,
            time_explosion=config.supernova.time_explosion,
            legacy_mode_enabled=legacy_mode_enabled,
        )
    else:
        packet_source = BlackBodySimpleSource(
            base_seed=config.montecarlo.seed,
            legacy_mode_enabled=legacy_mode_enabled,
        )

    return initialize_packet_source(packet_source, config, geometry)
