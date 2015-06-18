import logging

import networkx as nx
from tardis.plasma.exceptions import PlasmaMissingModule, NotInitializedModule

import tempfile


logger = logging.getLogger(__name__)


class BasePlasma(object):

    def __init__(self, plasma_modules, **kwargs):
        self.module_dict = {}
        self.input_modules = []
        self._init_modules(plasma_modules, **kwargs)
        self._build_graph()
        self.update(**kwargs)


    def __getattr__(self, item):
        if item in self.module_dict:
            return self.get_value(item)
        else:
            super(BasePlasma, self).__getattribute__(item)


    def __setattr__(self, key, value):
        if key != 'module_dict' and key in self.module_dict:
            raise AttributeError('Plasma inputs can only be updated using '
                                 'the \'update\' method')
        else:
            super(BasePlasma, self).__setattr__(key, value)

    def __dir__(self):
        attrs = [item for item in self.__dict__
                 if not item.startswith('_')]
        attrs += [item for item in self.__class__.__dict__
                 if not item.startswith('_')]
        attrs += self.module_dict.keys()

        return attrs

    def get_value(self, item):
        return self.module_dict[item].value

    def _build_graph(self):
        """
        Builds the directed Graph using network X

        :param plasma_modules:
        :return:
        """

        self.graph = nx.DiGraph()

        ## Adding all nodes
        self.graph.add_nodes_from([(key, {})
                                   for key, value in self.module_dict.items()])

        #Flagging all input modules
        self.input_modules = [key for key, item in self.module_dict.items()
                              if not hasattr(item, 'inputs')]

        for plasma_module in self.module_dict.values():
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

    def _update_module_type_str(self):
        for node in self.graph:
            self.module_dict[node]._update_type_str()

    def _resolve_update_list(self, changed_modules):
        """
        Returns a list of all plasma models which are affected by the
        changed_modules due to there dependency in the
        the plasma_graph.

        Parameters
        ----------

        graph: ~networkx.Graph
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

    def write_to_dot(self, fname):
        self._update_module_type_str()

        try:
            import pygraphviz
        except ImportError:
            raise ImportError('pygraphviz is needed for method '
                              '\'write_to_dot\'')

        for node in self.graph:
            self.graph.node[node]['label'] = self.module_dict[node].get_latex_label()
            self.graph.node[node]['color'] = 'red'
            self.graph.node[node]['shape'] = 'box '

        nx.write_dot(self.graph, fname)

    def write_to_tex(self, fname):
        try:
            import dot2tex
        except ImportError:
            raise ImportError('dot2tex is needed for method\'write_to_tex\'')

        temp_fname = tempfile.NamedTemporaryFile().name

        self.write_to_dot(temp_fname)

        dot_string = open(temp_fname).read()

        open(fname, 'w').write(dot2tex.dot2tex(dot_string, texmode='raw'))



class StandardPlasma(BasePlasma):

    def __init__(self, number_densities, atom_data, time_explosion,
                 delta_treatment=None, nlte_config=None, ionization_mode='lte',
                 excitation_mode='lte', w=None,
                 link_t_rad_t_electron=0.9):

        pass

