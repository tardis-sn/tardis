class Radial1DModel(object):
    def __init__(self, velocity, density, abundance, t_rad,
                 v_boundary_inner, v_boundary_outer, w=None):
        self.velocity = velocity
        self.density = density
        self.abundance = abundance
        self.t_rad = t_rad
        self.v_boundary_inner = v_boundary_inner
        self.v_boundary_outer = v_boundary_outer
        self.w = w

    @classmethod
    def from_config(cls, config):
        pass
