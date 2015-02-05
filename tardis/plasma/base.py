import networkx as nx
from tardis.plasma.exceptions import PlasmaMissingModule


class BasePlasma(object):

    input_modules = None
    
    def __init__(self, plasma_modules):

        self._build_dictionary(plasma_modules)
        self._build_graph(plasma_modules)

    def __getattr__(self, item):
        if item in self.input_modules:
            return self.input_modules[item].value
        else:
            super(BasePlasma).__getattribute__(item)

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
        self.input_modules = [item for item in plasma_modules
                              if not hasattr(item, 'inputs')]

        for plasma_module in plasma_modules:

            #Skipping any module that is an input module
            if plasma_module in self.input_modules:
                continue

            for input in plasma_module.inputs:
                if input not in self.graph:
                    raise PlasmaMissingModule('Module {0} requires input '
                                              '{1} which has not been added'
                                              ' to this plasma'.format(
                        plasma_module.name, input))
                self.graph.add_edge(input, plasma_module.name)

            



        def _build_dictionary(self, plasma_module):
            """
            Builds a dictionary with the plasma module names as keys
            :param plasma_module:
            :return:
            """
            self.module_dict = dict([
                (module.name, module) for module in plasma_module])


        def update_plasma(self, **kwargs):
            for key, value in kwargs.items():
                pass





class StandardPlasma(BasePlasma):

    def __init__(self, number_densities, atom_data, time_explosion,
                 delta_treatment=None, nlte_config=None, ionization_mode='lte',
                 excitation_mode='lte'):

        pass

