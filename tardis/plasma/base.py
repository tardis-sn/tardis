import os
import re
import logging
import tempfile
import fileinput

import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

from tardis.plasma.exceptions import PlasmaMissingModule, NotInitializedModule
from tardis.plasma.properties.base import *
from tardis.io.util import PlasmaWriterMixin

logger = logging.getLogger(__name__)

class BasePlasma(PlasmaWriterMixin):

    outputs_dict = {}
    hdf_name = 'plasma'
    def __init__(self, plasma_properties, property_kwargs=None, **kwargs):
        self.outputs_dict = {}
        self.input_properties = []
        self.plasma_properties = self._init_properties(plasma_properties,
                                                       property_kwargs, **kwargs)
        self._build_graph()
#        self.write_to_tex('Plasma_Graph')
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
        attrs += self.outputs_dict.keys()
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
                try:
                    position = self.outputs_dict[input].outputs.index(input)
                    label = self.outputs_dict[input].latex_name[position]
                    label = '$' + label + '$'
                except:
                    label = input.replace('_', '-')
                self.graph.add_edge(self.outputs_dict[input].name,
                    plasma_property.name, label = label)

    def _init_properties(self, plasma_properties, property_kwargs=None, **kwargs):
        """
        Builds a dictionary with the plasma module names as keys

        Parameters
        ----------

        plasma_modules: ~list
            list of Plasma properties
        property_kwargs: ~dict
            dict of plasma module : kwargs pairs. kwargs should be a dict
            of arguments that will be passed to the __init__ method of
            the respective plasma module.
        kwargs: dictionary
            input values for input properties. For example, t_rad=[5000, 6000,],
            j_blues=[..]

        """
        if property_kwargs is None:
            property_kwargs = {}
        plasma_property_objects = []
        self.previous_iteration_properties = []
        self.outputs_dict = {}
        for plasma_property in plasma_properties:

            if issubclass(plasma_property, PreviousIterationProperty):
                current_property_object = plasma_property(
                    **property_kwargs.get(plasma_property, {}))
                current_property_object.set_initial_value(kwargs)
                self.previous_iteration_properties.append(
                    current_property_object)

            elif issubclass(plasma_property, Input):
                if not set(kwargs.keys()).issuperset(plasma_property.outputs):
                    missing_input_values = (set(plasma_property.outputs) -
                    set(kwargs.keys()))
                    raise NotInitializedModule('Input {0} required for '
                                               'plasma but not given when '
                                               'instantiating the '
                                               'plasma'.format(
                                               missing_input_values))
                current_property_object = plasma_property(
                    **property_kwargs.get(plasma_property, {}))
            else:
                current_property_object = plasma_property(
                    self, **property_kwargs.get(plasma_property, {}))
            for output in plasma_property.outputs:
                self.outputs_dict[output] = current_property_object
                plasma_property_objects.append(current_property_object)
        return plasma_property_objects

    def store_previous_properties(self):
        for property in self.previous_iteration_properties:
            p = property.outputs[0]
            self.outputs_dict[p].set_value(
                self.get_value(re.sub(r'^previous_', '', p)))

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
        sort_order = list(nx.topological_sort(self.graph))

        descendants_ob.sort(key=lambda val: sort_order.index(val) )

        logger.debug('Updating modules in the following order: {}'.format(
            '->'.join(descendants_ob)))

        return descendants_ob

    def write_to_dot(self, fname, latex_label=True):
#        self._update_module_type_str()

        try:
            import pygraphviz
        except:
            logger.warn('pygraphviz missing. Plasma graph will not be '
                        'generated.')
            return
        print_graph = self.graph.copy()
        print_graph = self.remove_hidden_properties(print_graph)
        for node in print_graph:
            print_graph._node[str(node)]['label'] = node
            if hasattr(self.plasma_properties_dict[node],
                'latex_formula'):
                formulae = self.plasma_properties_dict[node].latex_formula
                for output in range(0, len(formulae)):
                    formula = formulae[output]
                    #label = formula.replace('\\', '\\\\')
                    print_graph._node[str(node)]['label']+='$'
                    print_graph._node[str(node)]['label']+=formula
                    print_graph._node[str(node)]['label']+='$'

        nx.drawing.nx_agraph.write_dot(print_graph, fname)

    def write_to_tex(self, fname_graph):
        try:
            import dot2tex
        except:
            logger.warn('dot2tex missing. Plasma graph will not be '
                        'generated.')
            return

        temp_fname = tempfile.NamedTemporaryFile().name

        self.write_to_dot(temp_fname)

        dot_string = open(temp_fname).read()

        open(fname_graph, 'w').write(dot2tex.dot2tex(dot_string,texmode='raw',format='tikz',tikzedgelabels=True))


        for line in fileinput.input(fname_graph, inplace = 1):
            print(line.replace(r'\documentclass{article}',
                r'\documentclass[class=minimal,border=20pt]{standalone}'),
                  end='')

        for line in fileinput.input(fname_graph, inplace = 1):
            print(line.replace(r'\enlargethispage{100cm}', ''), end='')

    def display_plasma_graph(self,node_size=2000,node_color='red',alpha=0.5,
                                font_size=7,edge_labels_font_color='blue',connectionstyle=None,node_shape='o',edge_font_size=4,pos='dot'):
        """
        TARDIS calculates the state of the gas of the exploding star in the plasma calculation. When the gas changes temperature, 
        TARDIS can automatically calculate what properties of the gas change. It uses a network graph for this. This also allows us 
        to visualize how the properties of the gas are calculated.
        This method displays an interactive DiGraph which lets us see this process.

        Parameters :
        node_size : Default - 2000, this value sets the size of node
        node_color : Default - 'red', Sets the color in which the node is displayed.
        alpha : Default - 0.5, sets the transperancy level of node
        font_size : Default - 7, font size of node and edge labels
        edge_labels_font_color : Default - 'blue', Color of edge labels
        edge_font_size : Default - 4, specify font size for edge labels
        connectionstyle : Default : None, Specifies the shape of connecting arrows. Refer https://matplotlib.org/3.1.1/gallery/userdemo/connectionstyle_demo.html for more details.

        Reference documentation : https://networkx.github.io/documentation/stable/reference/generated/networkx.drawing.nx_pylab.draw_networkx_labels.html#networkx.drawing.nx_pylab.draw_networkx_labels
        """

        try:
            import pygraphviz as pgv
        except:
            logger.warn('pygraphviz missing. Plasma graph will not be '
                        'generated.')

        temp_fname = tempfile.NamedTemporaryFile().name
        self.write_to_dot(temp_fname)
        G = nx.DiGraph(pgv.AGraph(temp_fname))

        nodes = G.nodes.data('label')._nodes
        mapping = {}
        for node in nodes:
            mapping[node] = nodes[node]['label']
        nx.relabel_nodes(G, mapping, copy=False)
        edge_labels = {}
        for u, v, d in G.edges(data=True):
            edge_labels[(u,v)] = d['label']

        try:
            layout = nx.nx_agraph.graphviz_layout(G, prog=pos)
        except :
            print("Unable to get layout {}, defaulting to 'dot'".format(pos))
            layout = nx.nx_agraph.graphviz_layout(G, prog='dot')

        fig = plt.gcf()
        fig.set_size_inches(25,15)
        edge_label_plot = nx.draw_networkx_edge_labels(G,layout,edge_labels=edge_labels,font_size=edge_font_size,font_color=edge_labels_font_color)
        nx.draw_networkx_nodes(G,layout,node_size=node_size,node_color=node_color,alpha=alpha, node_shape=node_shape, min_source_margin=3,min_target_margin=3)
        nx.draw_networkx_edges(G,layout,width=0.8,edge_color='black',connectionstyle=connectionstyle,arrows=True)
        nx.draw_networkx_labels(G,layout,font_size=font_size,font_family='sans-serif',)
        plt.axis('off')
        plt.show()

    def remove_hidden_properties(self, print_graph):
        for item in self.plasma_properties_dict.values():
            module = self.plasma_properties_dict[item.name].__class__
            if (issubclass(module, HiddenPlasmaProperty)):
                output = module.outputs[0]
                for value in self.plasma_properties_dict.keys():
                    if output in getattr(
                            self.plasma_properties_dict[value], 'inputs', []):
                        for input in self.plasma_properties_dict[
                            item.name].inputs:
                            try:
                                position = self.outputs_dict[
                                    input].outputs.index(input)
                                label = self.outputs_dict[
                                    input].latex_name[position]
                                label = '$' + label + '$'
                            except:
                                label = input.replace('_', '-')
                            self.graph.add_edge(self.outputs_dict[input].name,
                                value, label = label)
                print_graph.remove_node(str(item.name))
        return print_graph
