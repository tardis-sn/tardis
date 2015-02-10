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


class Abundance(DynamicInput):
    name = 'abundance'


class Density(DynamicInput):
    name = 'density'
    latex_name = r'$\rho$'

class TimeExplosion(DynamicInput):
    name = 'time_explosion'