from tardis.plasma.plasma_properties import BasePlasmaProperty



class PlasmaInput(object):

    def __init__(self):
        self.value = None
    def set_value(self, value):
        self.value = value

class StaticPlasmaInput(PlasmaInput):
    pass

class DynamicPlasmaInput(PlasmaInput):
    pass




class TRadiative(DynamicPlasmaInput):
    name = 't_rad'

class AtomicData(StaticPlasmaInput):
    name = 'atomic_data'
