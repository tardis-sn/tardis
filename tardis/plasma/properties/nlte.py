from tardis.plasma.properties.base import (BasePlasmaProperty,
                                           ProcessingPlasmaProperty)

__all__ = ['PreviousElectronDensities', 'PreviousBetaSobolevs',
           'HeliumNLTE']

class PreviousIterationProperty(BasePlasmaProperty):
    def _set_output_value(self, output, value):
        setattr(self, output, value)

    def set_value(self, value):
        assert len(self.outputs) == 1
        self._set_output_value(self.outputs[0], value)

class PreviousElectronDensities(PreviousIterationProperty):
    outputs = ('previous_electron_densities',)

class PreviousBetaSobolevs(PreviousIterationProperty):
    outputs = ('previous_beta_sobolevs',)

class HeliumNLTE(ProcessingPlasmaProperty):
    outputs = ('helium_population',)

    def calculate(self, level_boltzmann_factor, electron_densities,
        ionization_data, beta_rad, g, g_electron, w, t_rad, t_electrons,
        delta, zeta_data):
        pass
