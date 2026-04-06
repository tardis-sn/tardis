import numpy as np
from numba import njit, objmode, prange
from numba.np.ufunc.parallel import get_num_threads, get_thread_id
from numba.typed import List

from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.montecarlo import njit_dict
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
)
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    create_estimators_bulk_list,
    init_estimators_bulk,
)
from tardis.transport.montecarlo.estimators.estimators_continuum import (
    create_estimators_continuum_list,
    init_estimators_continuum,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    create_estimators_line_list,
    init_estimators_line,
)
from tardis.transport.montecarlo.modes.iip.packet_propagation import (
    packet_propagation,
)
from tardis.transport.montecarlo.packets.packet_collections import (
    PacketCollection,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    PacketStatus,
    RPacket,
)
from tardis.transport.montecarlo.progress_bars import update_packets_pbar


@njit(**njit_dict)
def montecarlo_transport(
    packet_collection: PacketCollection,
    geometry_state_numba: NumbaRadial1DGeometry,
    time_explosion: float,
    opacity_state_numba: OpacityStateNumba,
    montecarlo_configuration: MonteCarloConfiguration,
    n_levels_bf_species_by_n_cells_tuple: tuple,
    trackers: List,
    show_progress_bars: bool,
):
    """
    Main loop of the Monte Carlo radiative transfer routine for IIP mode.

    This function generates packet objects from the packet collection and
    propagates them through the ejecta, performing interactions with both
    lines and continuum processes, and collecting statistics for the
    radiative transfer simulation.

    Parameters
    ----------
    packet_collection : PacketCollection
        Collection containing initial packet properties (positions, directions,
        frequencies, energies, and seeds).
    geometry_state_numba : NumbaRadial1DGeometry
        Numba-compiled simulation geometry containing shell boundaries
        and velocity information.
    time_explosion : float
        Time since explosion in seconds, used for relativistic calculations.
    opacity_state_numba : OpacityStateNumba
        Numba-compiled opacity state containing line opacities, continuum
        opacities, and atomic data required for interactions.
    montecarlo_configuration : MonteCarloConfiguration
        Configuration object containing Monte Carlo simulation parameters
        and flags for various physics modules.
    n_levels_bf_species_by_n_cells_tuple : tuple
        Shape tuple for bound-free transitions (n_levels_bf_species, n_cells).
    trackers : numba.typed.List
        List of packet trackers for detailed packet interaction logging.
    show_progress_bars : bool
        Flag to enable/disable progress bar updates during simulation.

    Returns
    -------
    tuple
        A tuple containing:
        - estimators_bulk : EstimatorsBulk
            Updated bulk radiation field estimator object containing cell-level
            statistics collected during packet propagation.
        - estimators_line : EstimatorsLine
            Updated line radiation field estimator object containing line
            interaction statistics collected during packet propagation.
        - estimators_continuum : EstimatorsContinuum
            Updated continuum estimator object containing continuum interaction
            statistics collected during packet propagation.
    """
    no_of_packets = len(packet_collection.initial_nus)

    # Get the ID of the main thread and the number of threads
    main_thread_id = get_thread_id()
    n_threads = get_num_threads()

    # Note that get_thread_id() returns values from 0 to n_threads-1,
    # so we iterate from 0 to n_threads-1 to create the estimator_lists

    # Initialize estimators
    n_lines_by_n_cells_tuple = opacity_state_numba.tau_sobolev.shape
    n_cells = len(geometry_state_numba.r_inner)
    estimators_bulk = init_estimators_bulk(n_cells)
    estimators_line = init_estimators_line(n_lines_by_n_cells_tuple)
    estimators_continuum = init_estimators_continuum(
        n_levels_bf_species_by_n_cells_tuple, n_cells
    )

    # Initialize thread-local estimators
    estimators_bulk_list_thread = create_estimators_bulk_list(
        n_cells, n_threads
    )
    estimators_line_list_thread = create_estimators_line_list(
        n_lines_by_n_cells_tuple, n_threads
    )
    estimators_continuum_list_thread = create_estimators_continuum_list(
        n_levels_bf_species_by_n_cells_tuple, n_cells, n_threads
    )

    for i in prange(no_of_packets):
        thread_id = get_thread_id()
        if show_progress_bars:
            if thread_id == main_thread_id:
                with objmode:
                    update_amount = 1 * n_threads
                    update_packets_pbar(
                        update_amount,
                        no_of_packets,
                    )

        r_packet = RPacket(
            packet_collection.initial_radii[i],
            packet_collection.initial_mus[i],
            packet_collection.initial_nus[i],
            packet_collection.initial_energies[i],
            packet_collection.packet_seeds[i],
            i,
        )
        # Seed the random number generator
        np.random.seed(r_packet.seed)

        # Get the thread-local estimators for this thread
        estimators_bulk_thread = estimators_bulk_list_thread[thread_id]
        estimators_line_thread = estimators_line_list_thread[thread_id]
        estimators_continuum_thread = estimators_continuum_list_thread[
            thread_id
        ]

        # Get the RPacket tracker for this thread
        tracker = trackers[i]

        loop = packet_propagation(
            r_packet,
            geometry_state_numba,
            time_explosion,
            opacity_state_numba,
            estimators_bulk_thread,
            estimators_line_thread,
            estimators_continuum_thread,
            tracker,
            montecarlo_configuration,
        )
        packet_collection.output_nus[i] = r_packet.nu

        if r_packet.status == PacketStatus.REABSORBED:
            packet_collection.output_energies[i] = -r_packet.energy

        elif r_packet.status == PacketStatus.EMITTED:
            packet_collection.output_energies[i] = r_packet.energy

        # Finalize the tracker (e.g. trim arrays to actual size)
        tracker.finalize()

    for estimator_thread in estimators_bulk_list_thread:
        estimators_bulk.increment(estimator_thread)

    for estimator_thread in estimators_line_list_thread:
        estimators_line.increment(estimator_thread)

    for estimator_thread in estimators_continuum_list_thread:
        estimators_continuum.increment(estimator_thread)

    return (
        estimators_bulk,
        estimators_line,
        estimators_continuum,
    )
