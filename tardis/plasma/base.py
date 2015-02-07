import logging

import networkx as nx
from tardis.plasma.exceptions import PlasmaMissingModule, NotInitializedModule
from plasma_input import PlasmaInput

logger = logging.getLogger(__name__)


class BasePlasma(object):

    def __init__(self, plasma_modules, **kwargs):
        self.module_dict = {}
        self.input_modules = []
        self._init_modules(plasma_modules, **kwargs)
        self._build_graph(plasma_modules)
        self.update(**kwargs)

    def __getattr__(self, item):
        if item in self.module_dict:
            return self.module_dict[item].value
        else:
            super(BasePlasma, self).__getattribute__(item)


    def __setattr__(self, key, value):
        if key != 'module_dict' and key in self.module_dict:
            raise AttributeError('Plasma inputs can only be updated using '
                                 'the \'update\' method')
        else:
            super(BasePlasma, self).__setattr__(key, value)


    def _build_graph(self, plasma_modules):
        """
        Builds the directed Graph using network X

        :param plasma_modules:
        :return:
        """

        self.graph = nx.DiGraph()

        ## Adding all nodes
        self.graph.add_nodes_from([(key, {'label': value.get_label()})
                                   for key, value in self.module_dict.items()])

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
            if hasattr(module, 'set_value'):
                if module.name not in kwargs:
                    raise NotInitializedModule('Input {0} required for '
                                               'plasma but not given when '
                                               'instantiating the '
                                               'plasma'.format(module.name))
                current_module_object = module()
            else:
                current_module_object = module(self)

            self.module_dict[module.name] = current_module_object



    def update(self, **kwargs):
        for key in kwargs:
            if key not in self.module_dict:
                raise PlasmaMissingModule('Trying to update property {0}'
                                          ' that is unavailable'.format(key))
            self.module_dict[key].set_value(kwargs[key])


        for module_name in self._resolve_update_list(kwargs.keys()):

            self.module_dict[module_name].update()

    def _resolve_update_list(self, changed_modules):
        """
        Returns a list of all plasma models which are affected by the
        changed_modules due to there dependency in the
        the plasma_graph.

        Parameters
        ----------

        graph: ~networkx.Grapht
            the plasma graph as
        changed_modules: ~list
            all modules changed in the plasma

        Returns
        -------

            : ~list
            all affected modules.
        """

        descendants_ob = []

        for module in changed_modules:
            descendants_ob += nx.descendants(self.graph, module)

        descendants_ob = list(set(descendants_ob))
        sort_order = nx.topological_sort(self.graph)

        descendants_ob.sort(key=lambda val: sort_order.index(val) )

        logger.debug('Updating modules in the following order:'.format(
            '->'.join(descendants_ob)))

        return descendants_ob

    def write_to_dot(self):
        try:
            import pygraphviz
        except ImportError:
            raise ImportError('pygraphviz is needed for method '
                              '\'plasma_to_dot\'')




class StandardPlasma(BasePlasma):

    def __init__(self, number_densities, atom_data, time_explosion,
                 delta_treatment=None, nlte_config=None, ionization_mode='lte',
                 excitation_mode='lte'):

        pass

