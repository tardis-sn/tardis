from tardis.plasma.base_properties import BasePlasmaProperty



class Input(BasePlasmaProperty):

    def set_value(self, value):
        self.value = value


class StaticInput(Input):
    pass

class DynamicInput(Input):
    pass

class TRadiative(DynamicInput):
    name = 't_rad'
    latex_name = r'$T_\textrm{rad}$'


class DilutionFactor(DynamicInput):
    name = 'w'
    latex_name = r'$W$'


class AtomicData(StaticInput):
    name = 'atomic_data'
    latex_name = 'Atomic Data'

class Abundance(DynamicInput):
    name = 'abundance'
    latex_name = 'Abundance'
