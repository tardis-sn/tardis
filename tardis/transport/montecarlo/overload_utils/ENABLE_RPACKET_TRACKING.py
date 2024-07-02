import numba
from numba.typed import List

from tardis.transport.montecarlo.packet_trackers import RPacketTracker


def initialize_rpacket_trackers(
    ENABLE_RPACKET_TRACKING, no_of_packets, INITIAL_ARRAY_TRACKING_LENGTH
):
    pass


@numba.extending.overload(initialize_rpacket_trackers)
def ol_initialize_rpacket_trackers(
    ENABLE_RPACKET_TRACKING, no_of_packets, INITIAL_ARRAY_TRACKING_LENGTH
):
    if ENABLE_RPACKET_TRACKING.literal_value:

        def initialize(
            ENABLE_RPACKET_TRACKING,
            no_of_packets,
            INITIAL_ARRAY_TRACKING_LENGTH,
        ):
            rpacket_trackers = List()
            for i in range(no_of_packets):
                rpacket_trackers.append(
                    RPacketTracker(INITIAL_ARRAY_TRACKING_LENGTH)
                )
            return rpacket_trackers

        return initialize
    else:
        return (
            lambda ENABLE_RPACKET_TRACKING, no_of_packets, INITIAL_ARRAY_TRACKING_LENGTH: None
        )


def get_rpacket_tracker(ENABLE_RPACKET_TRACKING, rpacket_trackers, index):
    pass


@numba.extending.overload(get_rpacket_tracker)
def ol_get_rpacket_tracker(ENABLE_RPACKET_TRACKING, rpacket_trackers, index):
    if ENABLE_RPACKET_TRACKING.literal_value:
        return lambda ENABLE_RPACKET_TRACKING, rpacket_trackers, index: rpacket_trackers[
            index
        ]
    else:
        return lambda ENABLE_RPACKET_TRACKING, rpacket_trackers, index: None


def finalize_rpacket_trackers(ENABLE_RPACKET_TRACKING, rpacket_trackers):
    pass


@numba.extending.overload(finalize_rpacket_trackers)
def ol_finalize_rpacket_trackers(ENABLE_RPACKET_TRACKING, rpacket_trackers):
    if ENABLE_RPACKET_TRACKING.literal_value:

        def finalize(ENABLE_RPACKET_TRACKING, rpacket_trackers):
            for rpacket_tracker in rpacket_trackers:
                rpacket_tracker.finalize_array()
            return None

        return finalize
    else:
        return lambda ENABLE_RPACKET_TRACKING, rpacket_trackers: None


def track_rpacket(ENABLE_RPACKET_TRACKING, rpacket_tracker, r_packet):
    pass


@numba.extending.overload(track_rpacket)
def ol_track_rpacket(ENABLE_RPACKET_TRACKING, rpacket_tracker, r_packet):
    if ENABLE_RPACKET_TRACKING.literal_value:

        def track(ENABLE_RPACKET_TRACKING, rpacket_tracker, r_packet):
            rpacket_tracker.track(r_packet)
            return None

        return track
    else:
        return lambda ENABLE_RPACKET_TRACKING, rpacket_tracker, r_packet: None
