

class BasePlasma(object):
        pass


class PlasmaMissingModule(Exception):
    pass

class StandardPlasma(BasePlasma):

    def __init__(self, number_densities, atom_data, time_explosion,
                 delta_treatment=None, nlte_config=None, ionization_mode='lte',
                 excitation_mode='lte'):

        self.number_densities = number_densities
        self.atom_data = atom_data
        self.time_explosion = time_explosion
        self.nlte_config = nlte_config
        self.delta_treatment = delta_treatment
        self.electron_densities = self.number_densities.sum(axis=0)

    module_dict = {'partition_function':pf, 'level_populations':lp}

    def setup(self):
        for key in self.module_dict:
            current_module = self.module_dict[key]
            current_dependencies = current_module.get_dependencies()
            for key in current_dependencies:
                if key not in self.module_dict:
                    raise PlasmaMissingModule('Module {0} missing')
                current_module.register_depency()
    def evaluate(self):


