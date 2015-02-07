from tardis.plasma.plasma_properties import BasePlasmaProperty



class PlasmaInput(object):

    def __init__(self):
        self.value = None
    def set_value(self, value):
        self.value = value
    def get_label(self):
        return "Name: {0}\nType: {1}\n{2}".format(self.name, self.type_str,
                                                  getattr(self,
                                                          'latex_str', ''))



class StaticPlasmaInput(PlasmaInput):
    pass

class DynamicPlasmaInput(PlasmaInput):
    pass




class TRadiative(DynamicPlasmaInput):
    name = 't_rad'
    type_str = 'numpy.array'

class AtomicData(StaticPlasmaInput):
    name = 'atomic_data'
    type_str = 'tardis.atomic.AtomicData'
    label = 'Atomic Data'

class Abundance(DynamicPlasmaInput):
    name = 'abundance'
    type_str = 'pandas.DataFrame'