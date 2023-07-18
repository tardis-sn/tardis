from tardis.montecarlo.packet_source import BlackBodySimpleSource
from abc import ABC

class SourceFunction(ABC):
    pass

class PhotosphericBlackBody1D(SourceFunction):

    def __init__(self, v_inner_boundary, t_inner, no_of_packets, packet_source=None):

        self.packet_source = packet_source
        self.v_inner_boundary = v_inner_boundary
        self.t_inner = t_inner
        self.no_of_packets = no_of_packets