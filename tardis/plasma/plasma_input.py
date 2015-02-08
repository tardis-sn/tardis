from tardis.plasma.plasma_properties import BasePlasmaProperty



class Input(BasePlasmaProperty):

    def set_value(self, value):
        self.value = value

    def get_latex_label(self):
        return "Name: {0}".format(self.name)

class StaticInput(Input):
    pass

class DynamicInput(Input):
    pass

class TRadiative(DynamicInput):
    name = 't_rad'

class AtomicData(StaticInput):
    name = 'atomic_data'

class Abundance(DynamicInput):
    name = 'abundance'
