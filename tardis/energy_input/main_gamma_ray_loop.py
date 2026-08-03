import logging

import astropy.units as u
import numpy as np
import pandas as pd
from numba.typed import List as TypedList
from numpy.typing import NDArray  # noqa: TC002

from tardis.configuration.sorting_globals import SORTING_ALGORITHM
from tardis.energy_input.gamma_ray_transport import (
    calculate_ejecta_velocity_volume,
    iron_group_fraction_per_shell,
)
from tardis.energy_input.transport.gamma_packet_loop import gamma_packet_loop
from tardis.energy_input.transport.GXPacket import GXPacket
from tardis.energy_input.util import get_index
from tardis.io.atom_data import AtomData
from tardis.model.base import SimulationState
from tardis.transport.montecarlo.packet_source.high_energy import (
    GammaRayPacketSource,
    legacy_calculate_positron_fraction,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def calculate_electron_number_density(
    simulation_state: SimulationState,
    ejecta_volume: NDArray[np.float64],
    effective_time_array: NDArray[np.float64],
    legacy: bool = False,
    legacy_atom_data: AtomData | None = None,
) -> NDArray[np.float64]:
    """Calculate the time-dependent electron number density.

    Parameters
    ----------
    simulation_state : SimulationState
        State containing the ejecta geometry and composition.
    ejecta_volume : numpy.ndarray
        Shell volumes in cubic centimeters at the simulation-state time.
    effective_time_array : numpy.ndarray
        Effective times in seconds at which to evaluate the density.
    legacy : bool, optional
        If ``True``, calculate the elemental number density through the legacy
        simulation-state interface.
    legacy_atom_data : AtomData or None, optional
        Atomic data supplying elemental masses for the legacy calculation.
        Required when ``legacy`` is ``True``.

    Returns
    -------
    numpy.ndarray
        Electron number density in inverse cubic centimeters, indexed by shell
        and effective time.

    Raises
    ------
    ValueError
        If legacy mode is requested without ``legacy_atom_data``.
    """
    ejecta_velocity_volume = calculate_ejecta_velocity_volume(simulation_state)

    inv_volume_time = (
        1.0 / ejecta_velocity_volume[:, np.newaxis]
    ) / effective_time_array**3.0

    # Calculate the elemental number density
    if not legacy:
        elemental_number_density = (
            simulation_state.composition.isotopic_number_density.groupby(
                "atomic_number"
            ).sum()
        )
    else:
        if legacy_atom_data is None:
            raise ValueError(
                "legacy_atom_data must be provided when legacy=True"
            )
        elemental_number_density = (
            simulation_state.calculate_elemental_number_density(
                legacy_atom_data.atom_data.mass
            )
        )

    # Electron number density
    electron_number_density = elemental_number_density.mul(
        elemental_number_density.index,
        axis=0,
    ).sum()
    electron_number = np.array(electron_number_density * ejecta_volume)

    # Evolve electron number and mass density with time
    electron_number_density_time = (
        electron_number[:, np.newaxis] * inv_volume_time
    )

    return electron_number_density_time


def get_effective_time_array(
    time_start: float,
    time_end: float,
    time_space: str,
    time_steps: int,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Create time-bin boundaries and representative effective times.

    Parameters
    ----------
    time_start : float
        Start time in days.
    time_end : float
        End time in days.
    time_space : str
        Time-bin spacing, either ``"linear"`` or ``"log"``.
    time_steps : int
        Number of time bins.

    Returns
    -------
    times : numpy.ndarray
        Time-bin boundaries in days. The array has ``time_steps + 1`` entries.
    effective_time_array : numpy.ndarray
        Representative time of each bin in days. Logarithmic bins use the
        geometric mean and linear bins use the arithmetic mean.

    Raises
    ------
    AssertionError
        If ``time_start`` is not smaller than ``time_end``.
    """
    assert time_start < time_end, "time_start must be smaller than time_end!"
    if time_space == "log":
        times = np.geomspace(time_start, time_end, time_steps + 1)
        effective_time_array = np.sqrt(times[:-1] * times[1:])
    else:
        times = np.linspace(time_start, time_end, time_steps + 1)
        effective_time_array = 0.5 * (times[:-1] + times[1:])

    return times, effective_time_array


def run_gamma_ray_loop(
    simulation_state: SimulationState,
    legacy_isotope_decacy_df: pd.DataFrame,
    cumulative_decays_df: pd.DataFrame,
    number_of_packets: int,
    times: NDArray[np.float64],
    effective_time_array: NDArray[np.float64],
    seed: int,
    positronium_fraction: float,
    spectrum_bins: int,
    grey_opacity: float,
    photoabsorption_opacity: str = "tardis",
    pair_creation_opacity: str = "tardis",
    legacy: bool = False,
    legacy_atom_data: AtomData | None = None,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    """Propagate gamma-ray packets through homologously expanding ejecta.

    Parameters
    ----------
    simulation_state : SimulationState
        State containing the ejecta geometry, density, and composition.
    legacy_isotope_decacy_df : pandas.DataFrame
        Radioactive-decay transition data used to compute packet energies and
        isotope-specific positron fractions.
    cumulative_decays_df : pd.DataFrame
        Time-dependent radioactive-decay data from which packets are sampled.
    number_of_packets : int
        Number of Monte Carlo packets to propagate.
    times : numpy.ndarray
        Time-bin boundaries in days.
    effective_time_array : numpy.ndarray
        Representative time of each time bin in days.
    seed : int
        Seed for the random number generator.
    positronium_fraction : float
        Fraction of positrons that form positronium.
    spectrum_bins : int
        Number of logarithmically spaced escaping-spectrum energy bins.
    grey_opacity : float
        Grey opacity in square centimeters per gram. A negative value enables
        the detailed interaction opacities.
    photoabsorption_opacity : {"kasen", "tardis"}, optional
        Photoabsorption opacity prescription used when ``grey_opacity`` is
        negative.
    pair_creation_opacity : {"artis", "tardis"}, optional
        Pair-creation opacity prescription used when ``grey_opacity`` is
        negative.
    legacy : bool, optional
        Whether to use the legacy elemental-density and packet-energy
        calculations.
    legacy_atom_data : AtomData or None, optional
        Atomic data used by the legacy elemental-density calculation. Required
        when ``legacy`` is ``True``.

    Returns
    -------
    escape_energy : pandas.DataFrame
        Escaping spectral luminosity, indexed by energy in keV with time-bin
        columns in seconds.
    escape_energy_cosi : pandas.DataFrame
        Escaping photon rate per energy bin, indexed by energy in keV with
        time-bin columns in seconds.
    packets_df_escaped : pandas.DataFrame
        Final packet diagnostics, including status, frequencies, energies, and
        shell number.
    gamma_ray_deposited_energy : pandas.DataFrame
        Gamma-ray energy deposited in each shell and time bin, in ergs.
    total_deposited_energy : pandas.DataFrame
        Gamma-ray plus positron energy deposition rate in each shell and time
        bin, in ergs per second.
    positron_energy_df : pandas.DataFrame
        Positron energy deposited in each shell and time bin, in ergs.
    """
    np.random.seed(seed)
    times = times * u.d.to(u.s)
    effective_time_array = effective_time_array * u.d.to(u.s)
    inner_velocities = simulation_state.v_inner.to("cm/s").value
    outer_velocities = simulation_state.v_outer.to("cm/s").value
    ejecta_volume = simulation_state.volume.to("cm^3").value
    shell_masses = simulation_state.volume * simulation_state.density
    number_of_shells = len(shell_masses)
    # TODO: decaying upto times[0]. raw_isotope_abundance is possibly not the best name
    isotopic_mass_fraction = (
        simulation_state.composition.isotopic_mass_fraction.sort_values(
            by=["atomic_number", "mass_number"],
            ascending=False,
            kind=SORTING_ALGORITHM,
        )
    )

    dt_array = np.diff(times)

    # Calculate electron number density evolution
    electron_number_density_time = calculate_electron_number_density(
        simulation_state,
        ejecta_volume,
        effective_time_array,
        legacy=legacy,
        legacy_atom_data=legacy_atom_data,
    )

    # Calculate mass density evolution
    ejecta_velocity_volume = calculate_ejecta_velocity_volume(simulation_state)
    inv_volume_time = (
        1.0 / ejecta_velocity_volume[:, np.newaxis]
    ) / effective_time_array**3.0
    mass_density_time = shell_masses[:, np.newaxis] * inv_volume_time

    packet_source = GammaRayPacketSource(
        cumulative_decays_df=cumulative_decays_df,
        isotope_decay_df=legacy_isotope_decacy_df,
        positronium_fraction=positronium_fraction,
        inner_velocities=inner_velocities,
        outer_velocities=outer_velocities,
        times=times,
        effective_times=effective_time_array,
        base_seed=seed,
    )

    logger.info("Creating packets")
    if legacy:
        # Calculate energy per packet for legacy mode using legacy_isotope_decacy_df
        gamma_df = legacy_isotope_decacy_df[
            legacy_isotope_decacy_df["radiation"] == "g"
        ]
        total_energy_gamma = gamma_df["decay_energy_erg"].sum()
        energy_per_packet = total_energy_gamma / number_of_packets
        legacy_energy_per_packet = energy_per_packet
    else:
        # Let the packet source calculate energy per packet internally
        legacy_energy_per_packet = None
        energy_per_packet = None

    packet_collection = packet_source.create_packets(
        cumulative_decays_df,
        number_of_packets,
        legacy_energy_per_packet=legacy_energy_per_packet,
    )

    # For non-legacy mode, get the energy per packet from the packet source calculation
    if not legacy:
        gamma_df = cumulative_decays_df[
            cumulative_decays_df["radiation"] == "g"
        ]
        total_energy_gamma = gamma_df["decay_energy_erg"].sum()
        energy_per_packet = total_energy_gamma / number_of_packets

    total_energy = np.zeros((number_of_shells, len(times) - 1))

    logger.info("Creating packet list")
    # This for loop is expensive. Need to rewrite GX packet to handle arrays
    packets = TypedList()
    for i in range(number_of_packets):
        packets.append(
            GXPacket(
                packet_collection.location[:, i],
                packet_collection.direction[:, i],
                packet_collection.energy_rf[i],
                packet_collection.energy_cmf[i],
                packet_collection.nu_rf[i],
                packet_collection.nu_cmf[i],
                packet_collection.status[i],
                packet_collection.shell[i],
                packet_collection.time_start[i],
                packet_collection.time_index[i],
            )
        )

    # Calculate isotope positron fraction separately
    isotope_positron_fraction = legacy_calculate_positron_fraction(
        legacy_isotope_decacy_df,
        packet_collection.source_isotopes,
        number_of_packets,
    )
    for i, p in enumerate(packets):
        total_energy[p.shell, p.time_idx] += (
            isotope_positron_fraction[i] * energy_per_packet
        )

    logger.info(
        "Total energy deposited by the positrons is %s",
        total_energy.sum().sum(),
    )

    # positron_energy = total_energy

    # Copy the positron energy to a new dataframe
    positron_energy_df = pd.DataFrame(
        data=total_energy.copy(), columns=times[:-1]
    )

    energy_bins = np.logspace(2, 3.8, spectrum_bins)
    energy_out = np.zeros((len(energy_bins - 1), len(times) - 1))
    energy_out_cosi = np.zeros((len(energy_bins - 1), len(times) - 1))
    energy_deposited = np.zeros((number_of_shells, len(times) - 1))
    packets_info_array = np.zeros((int(number_of_packets), 8))
    iron_group_fraction = iron_group_fraction_per_shell(simulation_state)

    logger.info("Entering the main gamma-ray loop")

    total_cmf_energy = 0
    total_rf_energy = 0

    for p in packets:
        total_cmf_energy += p.energy_cmf
        total_rf_energy += p.energy_rf

    logger.info("Total CMF energy is %s", total_cmf_energy)
    logger.info("Total RF energy is %s", total_rf_energy)

    (
        energy_out,
        energy_out_cosi,
        packets_array,
        energy_deposited_gamma,
        total_energy,
    ) = gamma_packet_loop(
        packets,
        grey_opacity,
        photoabsorption_opacity,
        pair_creation_opacity,
        electron_number_density_time,
        mass_density_time,
        iron_group_fraction.to_numpy(),
        inner_velocities,
        outer_velocities,
        dt_array,
        times,
        effective_time_array,
        energy_bins,
        energy_out,
        energy_out_cosi,
        total_energy,
        energy_deposited,
        packets_info_array,
    )

    packets_df_escaped = pd.DataFrame(
        data=packets_array,
        columns=[
            "packet_index",
            "status",
            "nu_cmf",
            "nu_rf",
            "energy_cmf",
            "luminosity",
            "energy_rf",
            "shell_number",
        ],
    )

    escape_energy = pd.DataFrame(
        data=energy_out, columns=times[:-1], index=energy_bins
    )

    escape_energy_cosi = pd.DataFrame(
        data=energy_out_cosi, columns=times[:-1], index=energy_bins
    )

    # deposited energy by gamma-rays in ergs
    gamma_ray_deposited_energy = pd.DataFrame(
        data=energy_deposited_gamma, columns=times[:-1]
    )
    # deposited energy by positrons and gamma-rays in ergs
    total_energy = pd.DataFrame(data=total_energy, columns=times[:-1])

    logger.info(
        "Total energy deposited by gamma rays and positrons is %s",
        total_energy.sum().sum(),
    )

    total_deposited_energy = total_energy / dt_array

    return (
        escape_energy,
        escape_energy_cosi,
        packets_df_escaped,
        gamma_ray_deposited_energy,
        total_deposited_energy,
        positron_energy_df,
    )


def get_packet_properties(
    number_of_shells: int,
    times: NDArray[np.float64],
    time_steps: int,
    packets: list[GXPacket],
) -> tuple[
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
]:
    """Bin packet frequencies and energies by shell and time.

    Parameters
    ----------
    number_of_shells : int
        Number of ejecta shells.
    times : numpy.ndarray
        Time-bin boundaries in the same units as each packet's current time.
    time_steps : int
        Number of time bins.
    packets : list of GXPacket
        Gamma-ray packets to bin by shell and time.

    Returns
    -------
    packets_nu_cmf_array : numpy.ndarray
        Sum of comoving-frame frequencies in each shell and time bin.
    packets_nu_rf_array : numpy.ndarray
        Sum of rest-frame frequencies in each shell and time bin.
    packets_energy_cmf_array : numpy.ndarray
        Sum of comoving-frame energies in each shell and time bin.
    packets_energy_rf_array : numpy.ndarray
        Sum of rest-frame energies in each shell and time bin.
    packets_positron_energy_array : numpy.ndarray
        Sum of positron energies in each shell and time bin.
    """
    shell_number = []
    time_current = []

    # Bin the frequency of the packets in shell and time

    packets_nu_cmf_array = np.zeros((number_of_shells, time_steps))
    packets_nu_rf_array = np.zeros((number_of_shells, time_steps))
    packets_energy_cmf_array = np.zeros((number_of_shells, time_steps))
    packets_energy_rf_array = np.zeros((number_of_shells, time_steps))
    packets_positron_energy_array = np.zeros((number_of_shells, time_steps))

    for p in packets:
        time_index = get_index(p.time_current, times)
        shell_number.append(p.shell)
        time_current.append(p.time_current)
        packets_nu_cmf_array[p.shell, time_index] += p.nu_cmf
        packets_nu_rf_array[p.shell, time_index] += p.nu_rf
        packets_energy_cmf_array[p.shell, time_index] += p.energy_cmf
        packets_energy_rf_array[p.shell, time_index] += p.energy_rf
        packets_positron_energy_array[p.shell, time_index] += p.positron_energy

    return (
        packets_nu_cmf_array,
        packets_nu_rf_array,
        packets_energy_cmf_array,
        packets_energy_rf_array,
        packets_positron_energy_array,
    )
