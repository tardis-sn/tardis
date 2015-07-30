from tardis.plasma.properties.base import BasePlasmaProperty

__all__ = ['PreviousElectronDensities', 'PreviousBetaSobolevs']

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