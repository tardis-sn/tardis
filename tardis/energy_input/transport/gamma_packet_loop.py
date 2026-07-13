import numpy as np
from numba import njit
from numpy.typing import NDArray  # noqa: TC002

from tardis.energy_input.transport.gamma_ray_grid import (
    distance_trace,
    move_packet,
)
from tardis.energy_input.transport.gamma_ray_interactions import (
    compton_scatter,
    get_compton_fraction_artis,
    pair_creation_packet,
    scatter_type,
)
from tardis.energy_input.transport.GXPacket import GXPacket, GXPacketStatus
from tardis.energy_input.util import (
    C_CGS,
    H_CGS_KEV,
    doppler_factor_3d,
    get_index,
)
from tardis.opacities.opacities import (
    SIGMA_T,
    compton_opacity_calculation,
    kappa_calculation,
    pair_creation_opacity_artis,
    pair_creation_opacity_calculation,
    photoabsorption_opacity_calculation,
)
from tardis.transport.montecarlo import njit_dict_no_parallel


@njit(**njit_dict_no_parallel)
def gamma_packet_loop(
    packets: list[GXPacket],
    grey_opacity: float,
    photoabsorption_opacity_type: str,
    pair_creation_opacity_type: str,
    electron_number_density_time: NDArray[np.float64],
    mass_density_time: NDArray[np.float64],
    iron_group_fraction_per_shell: NDArray[np.float64],
    inner_velocities: NDArray[np.float64],
    outer_velocities: NDArray[np.float64],
    dt_array: NDArray[np.float64],
    times: NDArray[np.float64],
    effective_time_array: NDArray[np.float64],
    energy_bins: NDArray[np.float64],
    energy_out: NDArray[np.float64],
    energy_out_cosi: NDArray[np.float64],
    total_energy: NDArray[np.float64],
    energy_deposited_gamma: NDArray[np.float64],
    packets_info_array: NDArray[np.float64],
) -> tuple[
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
]:
    """Propagate gamma-ray packets through the simulation.

    Parameters
    ----------
    packets : list of GXPacket
        Packets to propagate through the ejecta.
    grey_opacity : float
        Grey opacity value in cm^2/g. Negative values trigger detailed opacity
        calculations.
    photoabsorption_opacity_type : str
        Method used to calculate photoabsorption opacity.
    pair_creation_opacity_type : str
        Method used to calculate pair creation opacity.
    electron_number_density_time : array float64
        Electron number densities indexed by shell and time.
    mass_density_time : array float64
        Mass densities indexed by shell and time.
    iron_group_fraction_per_shell : array float64
        Iron-group mass fraction for each shell.
    inner_velocities : array float64
        Inner shell velocities.
    outer_velocities : array float64
        Outer shell velocities.
    dt_array : array float64
        Simulation delta-time steps.
    times : array float64
        Simulation time-step boundaries.
    effective_time_array : array float64
        Effective simulation times for opacity and geometry calculations.
    energy_bins : array float64
        Energy bins for escaping gamma-rays.
    energy_out : array float64
        Escaped energy array updated in-place.
    energy_out_cosi : array float64
        Escaped photon count-rate array updated in-place.
    total_energy : array float64
        Total deposited energy array updated in-place.
    energy_deposited_gamma : array float64
        Gamma-ray deposited energy array updated in-place.
    packets_info_array : array float64
        Packet diagnostic output array updated in-place.

    Returns
    -------
    energy_out : array float64
        Escaped energy array.
    energy_out_cosi : array float64
        Escaped photon count-rate array.
    packets_info_array : array float64
        Packet diagnostic output array.
    energy_deposited_gamma : array float64
        Gamma-ray deposited energy array.
    total_energy : array float64
        Total deposited energy array.

    Raises
    ------
    ValueError
        Packet time index less than zero
    """
    escaped_packets = 0
    scattered_packets = 0
    packet_count = len(packets)
    # Logging does not work with numba. Using print instead.
    print("Entering gamma ray loop for " + str(packet_count) + " packets")

    for packet_idx in range(packet_count):
        packet = packets[packet_idx]
        time_idx = packet.time_idx

        if time_idx < 0:
            print(packet.time_start, time_idx)
            raise ValueError("Packet time index less than 0!")

        scattered = False
        # Not used now. Useful for the deposition estimator.
        # initial_energy = packet.energy_cmf

        while packet.status == GXPacketStatus.IN_PROCESS:
            # Get delta-time value for this step
            dt = dt_array[time_idx]
            # Calculate packet comoving energy for opacities
            comoving_energy = H_CGS_KEV * packet.nu_cmf

            if grey_opacity < 0:
                doppler_factor = doppler_factor_3d(
                    packet.direction,
                    packet.location,
                    times[time_idx],
                )

                kappa = kappa_calculation(comoving_energy)

                # artis threshold for Thomson scattering
                if kappa < 1e-2:
                    compton_opacity = (
                        SIGMA_T * electron_number_density_time[packet.shell, time_idx]
                    )
                else:
                    compton_opacity = compton_opacity_calculation(
                        comoving_energy,
                        electron_number_density_time[packet.shell, time_idx],
                    )

                if photoabsorption_opacity_type == "kasen":
                    # currently not functional, requires proton count and
                    # electron count per isotope
                    photoabsorption_opacity = 0
                    # photoabsorption_opacity_calculation_kasen()
                elif photoabsorption_opacity_type == "tardis":
                    photoabsorption_opacity = photoabsorption_opacity_calculation(
                        comoving_energy,
                        mass_density_time[packet.shell, time_idx],
                        iron_group_fraction_per_shell[packet.shell],
                    )
                else:
                    raise ValueError("Invalid photoabsorption opacity type!")

                if pair_creation_opacity_type == "artis":
                    pair_creation_opacity = pair_creation_opacity_artis(
                        comoving_energy,
                        mass_density_time[packet.shell, time_idx],
                        iron_group_fraction_per_shell[packet.shell],
                    )
                elif pair_creation_opacity_type == "tardis":
                    pair_creation_opacity = pair_creation_opacity_calculation(
                        comoving_energy,
                        mass_density_time[packet.shell, time_idx],
                        iron_group_fraction_per_shell[packet.shell],
                    )
                else:
                    raise ValueError("Invalid pair creation opacity type!")
            else:
                compton_opacity = 0.0
                pair_creation_opacity = 0.0
                photoabsorption_opacity = (
                    grey_opacity * mass_density_time[packet.shell, time_idx]
                )

            # convert opacities to rest frame
            total_opacity = (
                compton_opacity + photoabsorption_opacity + pair_creation_opacity
            ) * doppler_factor

            packet.tau = -np.log(np.random.random())

            (
                distance_interaction,
                distance_boundary,
                distance_time,
                shell_change,
            ) = distance_trace(
                packet,
                inner_velocities,
                outer_velocities,
                total_opacity,
                effective_time_array[time_idx],
                times[time_idx + 1],
            )

            distance = min(distance_interaction, distance_boundary, distance_time)

            packet.time_start += distance / C_CGS

            packet = move_packet(packet, distance)

            if distance == distance_time:
                time_idx += 1

                if time_idx > len(effective_time_array) - 1:
                    # Packet ran out of time
                    packet.status = GXPacketStatus.END
                else:
                    packet.shell = get_index(
                        packet.get_location_r(),
                        inner_velocities * times[time_idx],
                    )

            elif distance == distance_interaction:
                packet.status = scatter_type(
                    compton_opacity,
                    photoabsorption_opacity,
                    total_opacity,
                )

                packet, ejecta_energy_gained = process_packet_path(packet)

                # Ejecta gains energy from the packets (gamma-rays)
                energy_deposited_gamma[packet.shell, time_idx] += ejecta_energy_gained
                # Ejecta gains energy from both gamma-rays and positrons
                total_energy[packet.shell, time_idx] += ejecta_energy_gained

                if packet.status == GXPacketStatus.PHOTOABSORPTION:
                    # Packet destroyed, go to the next packet
                    break
                packet.status = GXPacketStatus.IN_PROCESS
                scattered = True

            else:
                packet.shell += shell_change

                if packet.shell > len(mass_density_time[:, 0]) - 1:
                    rest_energy = packet.nu_rf * H_CGS_KEV
                    energy_bin_idx = get_index(rest_energy, energy_bins)
                    bin_width = energy_bins[energy_bin_idx + 1] - energy_bins[energy_bin_idx]
                    freq_bin_width = bin_width / H_CGS_KEV

                    # get energy out in ergs per second per keV
                    energy_out[energy_bin_idx, time_idx] += (
                        packet.energy_rf
                        / dt
                        / freq_bin_width  # Take light crossing time into account
                    )
                    # get energy out in photons per second per keV
                    energy_out_cosi[energy_bin_idx, time_idx] += 1 / dt / bin_width

                    luminosity = packet.energy_rf / dt
                    packet.status = GXPacketStatus.ESCAPED
                    escaped_packets += 1
                    if scattered:
                        scattered_packets += 1
                elif packet.shell < 0:
                    packet.energy_rf = 0.0
                    packet.energy_cmf = 0.0
                    packet.status = GXPacketStatus.END

            packets_info_array[packet_idx] = np.array(
                [
                    packet_idx,
                    packet.status,
                    packet.nu_cmf,
                    packet.nu_rf,
                    packet.energy_cmf,
                    luminosity,
                    packet.energy_rf,
                    packet.shell,
                ]
            )

    print("Number of escaped packets:", escaped_packets)
    print("Number of scattered packets:", scattered_packets)

    return (
        energy_out,
        energy_out_cosi,
        packets_info_array,
        energy_deposited_gamma,
        total_energy,
    )


@njit(**njit_dict_no_parallel)
def process_packet_path(packet: GXPacket) -> tuple[GXPacket, float]:
    """Move a packet through interactions.

    Parameters
    ----------
    packet : GXPacket
        Packet to process.

    Returns
    -------
    packet : GXPacket
        Packet after processing.
    ejecta_energy_gained : float
        Energy injected into the ejecta.
    """
    if packet.status == GXPacketStatus.COMPTON_SCATTER:
        comoving_freq_energy = packet.nu_cmf * H_CGS_KEV

        compton_angle, compton_fraction = get_compton_fraction_artis(
            comoving_freq_energy
        )

        # Packet is no longer a gamma-ray, destroy it
        if np.random.random() < 1 / compton_fraction:
            packet.nu_cmf = packet.nu_cmf / compton_fraction

            packet.direction = compton_scatter(packet, compton_angle)

            # Calculate rest frame frequency after scaling by the fraction that remains
            doppler_factor = doppler_factor_3d(
                packet.direction,
                packet.location,
                packet.time_start,
            )

            packet.nu_rf = packet.nu_cmf / doppler_factor
            packet.energy_rf = packet.energy_cmf / doppler_factor

            ejecta_energy_gained = 0.0
        else:
            packet.status = GXPacketStatus.PHOTOABSORPTION

    if packet.status == GXPacketStatus.PAIR_CREATION:
        packet = pair_creation_packet(packet)
        ejecta_energy_gained = 0.0

    if packet.status == GXPacketStatus.PHOTOABSORPTION:
        # Ejecta gains comoving energy
        ejecta_energy_gained = packet.energy_cmf

    return packet, ejecta_energy_gained
