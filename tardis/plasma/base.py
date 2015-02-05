import networkx as nx
from tardis.plasma.exceptions import PlasmaMissingModule, NotInitializedModule
from plasma_input import PlasmaInput


class BasePlasma(object):

    input_modules = []

    def __init__(self, plasma_modules, **kwargs):

        self._init_modules(plasma_modules, **kwargs)
        self._build_graph(plasma_modules)

    def __getattr__(self, item):
        if item in self.module_dict:
            return self.module_dict[item].value
        else:
            super(BasePlasma, self).__getattribute__(item)

    def _build_graph(self, plasma_modules):
        """
        Builds the directed Graph using network X

        :param plasma_modules:
        :return:
        """

        self.graph = nx.DiGraph()

        ## Adding all nodes
        self.graph.add_nodes_from(self.module_dict.keys())

        #Flagging all input modules
        self.input_modules = [item.name for item in plasma_modules
                              if not hasattr(item, 'inputs')]

        for plasma_module in plasma_modules:

            #Skipping any module that is an input module
            if plasma_module.name in self.input_modules:
                continue

            for input in plasma_module.inputs:
                if input not in self.graph:
                    raise PlasmaMissingModule('Module {0} requires input '
                                              '{1} which has not been added'
                                              ' to this plasma'.format(
                        plasma_module.name, input))
                self.graph.add_edge(input, plasma_module.name)

    def _init_modules(self, plasma_modules, **kwargs):
        """
        Builds a dictionary with the plasma module names as keys
        :param plasma_modules:
        :return:
        """

        self.module_dict = {}
        for module in plasma_modules:
            if not hasattr(module, 'inputs'):
                if module.name not in kwargs:
                    raise NotInitializedModule('Input {0} required for '
                                               'plasma but not given when '
                                               'instantiating the '
                                               'plasma'.format(module.name))
                current_module_object = module(kwargs[module.name])
            else:
                current_module_object = module(self)

            self.module_dict[module.name] = current_module_object



    def update(self, **kwargs):
        for key in kwargs:
            if key not in self.module_dict:
                raise PlasmaMissingModule('Trying to update property {0}'
                                          ' that is unavailable'.format(key))

        update_list = self._resolve_update_list(kwargs.keys())

        for module in update_list:
            module.update()





class StandardPlasma(BasePlasma):

    def __init__(self, number_densities, atom_data, time_explosion,
                 delta_treatment=None, nlte_config=None, ionization_mode='lte',
                 excitation_mode='lte'):

        pass

