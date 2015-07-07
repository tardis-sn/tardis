import logging

import networkx as nx
from tardis.plasma.exceptions import PlasmaMissingModule, NotInitializedModule

import tempfile


logger = logging.getLogger(__name__)


class BasePlasma(object):
    outputs_dict = {}
    def __init__(self, plasma_properties, **kwargs):
        self.outputs_dict = {}
        self.input_properties = []
        self.plasma_properties = self._init_properties(plasma_properties,
                                                       **kwargs)
        
        self._build_graph()
        self.update(**kwargs)

    def __getattr__(self, item):
        if item in self.outputs_dict:
            return self.get_value(item)
        else:
            super(BasePlasma, self).__getattribute__(item)

    def __setattr__(self, key, value):
        if key != 'module_dict' and key in self.outputs_dict:
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

    @property
    def plasma_properties_dict(self):
        return {item.name:item for item in self.plasma_properties}

    def get_value(self, item):
        return getattr(self.outputs_dict[item], item)

    def _build_graph(self):
        """
        Builds the directed Graph using network X

        :param plasma_modules:
        :return:
        """

        self.graph = nx.DiGraph()

        ## Adding all nodes
        self.graph.add_nodes_from([(plasma_property.name, {})
                                   for plasma_property
                                   in self.plasma_properties])

        #Flagging all input modules
        self.input_properties = [item for item in self.plasma_properties
                                if not hasattr(item, 'inputs')]

        for plasma_property in self.plasma_properties:
            #Skipping any module that is an input module
            if plasma_property in self.input_properties:
                continue

            for input in plasma_property.inputs:
                if input not in self.outputs_dict:
                    raise PlasmaMissingModule('Module {0} requires input '
                                              '{1} which has not been added'
                                              ' to this plasma'.format(
                        plasma_property.name, input))
                self.graph.add_edge(self.outputs_dict[input].name,
                plasma_property.name, label=input)

    def _init_properties(self, plasma_properties, **kwargs):
        """
        Builds a dictionary with the plasma module names as keys

        Parameters
        ----------

        plasma_modules: ~list
            list of Plasma properties
        kwargs: dictionary
            input values for input properties. For example, t_rad=[5000, 6000,],
            j_blues=[..]

        """
        plasma_property_objects = []
        self.outputs_dict = {}
        for plasma_property in plasma_properties:

            if hasattr(plasma_property, 'set_value'):
                #duck-typing for PlasmaInputProperty
                #that means if it is an input property from model
                if not set(kwargs.keys()).issuperset(plasma_property.outputs):
                    missing_input_values = (set(plasma_property.outputs) -
                    set(kwargs.keys()))
                    raise NotInitializedModule('Input {0} required for '
                                               'plasma but not given when '
                                               'instantiating the '
                                               'plasma'.format(
                                               missing_input_values))
                current_property_object = plasma_property()
            else:
                current_property_object = plasma_property(self)
            for output in plasma_property.outputs:
                self.outputs_dict[output] = current_property_object
                plasma_property_objects.append(current_property_object)
        return plasma_property_objects

    def update(self, **kwargs):
        for key in kwargs:
            if key not in self.outputs_dict:
                raise PlasmaMissingModule('Trying to update property {0}'
                                          ' that is unavailable'.format(key))
            self.outputs_dict[key].set_value(kwargs[key])

        for module_name in self._resolve_update_list(kwargs.keys()):
            self.plasma_properties_dict[module_name].update()

    def _update_module_type_str(self):
        for node in self.graph:
            self.outputs_dict[node]._update_type_str()

    def _resolve_update_list(self, changed_properties):
        """
        Returns a list of all plasma models which are affected by the
        changed_modules due to there dependency in the
        the plasma_graph.

        Parameters
        ----------

        changed_modules: ~list
            all modules changed in the plasma

        Returns
        -------

            : ~list
            all affected modules.
        """

        descendants_ob = []

        for plasma_property in changed_properties:
            node_name = self.outputs_dict[plasma_property].name
            descendants_ob += nx.descendants(self.graph, node_name)

        descendants_ob = list(set(descendants_ob))
        sort_order = nx.topological_sort(self.graph)

        descendants_ob.sort(key=lambda val: sort_order.index(val) )

        logger.debug('Updating modules in the following order:'.format(
            '->'.join(descendants_ob)))

        return descendants_ob

    def write_to_dot(self, fname, latex_label=True):
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
                 link_t_rad_t_electron=0.9, nlte_species=None,
                 previous_beta_sobolevs=None):

        pass
