class PlasmaInput(object):
    pass


class TRadiative(PlasmaInput):

    name = 't_rad'

    def __init__(self, t_rad):
        self.t_rad = t_rad

    def calculate(self):
        return self.t_rad
