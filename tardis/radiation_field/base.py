import numpy as np

class RadiationField:
    """_summary_

    Parameters
    ----------
    t_rad : numpy.ndarray
        Radiative temperature in each shell 
    w : numpy.ndarray
        Dilution Factors in each shell
    opacities : Opacites
        Opacity container object
    source_function : SourceFunction
        Source function for radiative transfer, for example a packet_source
    """

    def __init__(self, t_rad, w, opacities, source_function):

        self.t_rad = t_rad
        self.w = w
        self.opacities = opacities
        self.source_function = source_function