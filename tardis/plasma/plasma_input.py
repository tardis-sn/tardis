from tardis.plasma.plasma_properties import BasePlasmaProperty



class PlasmaInput(BasePlasmaProperty):
    pass

class StaticPlasmaInput(PlasmaInput):
    pass

class DynamicPlasmaInput(PlasmaInput):
    pass




class TRadiative(DynamicPlasmaInput):

    name = 't_rad'

    def __init__(self, t_rad):
        self.value = t_rad


class AtomicData(StaticPlasmaInput):

    name = 'atomic_data'

    def __init__(self, atomic_data):
        self.value = atomic_data
