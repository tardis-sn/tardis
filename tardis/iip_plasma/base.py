import fileinput
import logging
import tempfile

import networkx as nx

from tardis.iip_plasma.exceptions import (
    NotInitializedModule,
    PlasmaMissingModule,
)
from tardis.iip_plasma.properties.base import *

logger = logging.getLogger(__name__)


class BasePlasma:
    outputs_dict = {}

    def __init__(self, plasma_properties, **kwargs):
        self.outputs_dict = {}
        self.input_properties = []
        self.plasma_properties = self._init_properties(
            plasma_properties, **kwargs
        )
        self._build_graph()
        #        self.write_to_tex('Plasma_Graph')
        self.plasma_converged = True
        self.update(**kwargs)

    def __getattr__(self, item):
        if item in self.outputs_dict:
            return self.get_value(item)
        super(BasePlasma, self).__getattribute__(item)

    def __setattr__(self, key, value):
        if key != "module_dict" and key in self.outputs_dict:
            raise AttributeError(
                "Plasma inputs can only be updated using the 'update' method"
            )
        super(BasePlasma, self).__setattr__(key, value)

    def __dir__(self):
        attrs = [item for item in self.__dict__ if not item.startswith("_")]
        attrs += [
            item for item in self.__class__.__dict__ if not item.startswith("_")
        ]
        attrs += self.module_dict.keys()
        return attrs

    @property
    def plasma_properties_dict(self):
        return {item.name: item for item in self.plasma_properties}

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
        self.graph.add_nodes_from(
            [
                (plasma_property.name, {})
                for plasma_property in self.plasma_properties
            ]
        )

        # Flagging all input modules
        self.input_properties = [
            item
            for item in self.plasma_properties
            if not hasattr(item, "inputs")
        ]

        for plasma_property in self.plasma_properties:
            # Skipping any module that is an input module
            if plasma_property in self.input_properties:
                continue

            for input in plasma_property.inputs:
                if input not in self.outputs_dict:
                    raise PlasmaMissingModule(
                        f"Module {plasma_property.name} requires input "
                        f"{input} which has not been added"
                        " to this plasma"
                    )
                try:
                    position = self.outputs_dict[input].outputs.index(input)
                    label = self.outputs_dict[input].latex_name[position]
                    label = "$" + label + "$"
                    label = label.replace("\\", "\\\\")
                except:
                    label = input.replace("_", "-")
                self.graph.add_edge(
                    self.outputs_dict[input].name,
                    plasma_property.name,
                    label=label,
                )

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
        self.previous_iteration_properties = []
        self.outputs_dict = {}
        for plasma_property in plasma_properties:
            if issubclass(plasma_property, PreviousIterationProperty):
                current_property_object = plasma_property()
                current_property_object.set_initial_value(kwargs)
                self.previous_iteration_properties.append(
                    current_property_object
                )

            elif issubclass(plasma_property, Input):
                if not set(kwargs.keys()).issuperset(plasma_property.outputs):
                    missing_input_values = set(plasma_property.outputs) - set(
                        kwargs.keys()
                    )
                    raise NotInitializedModule(
                        f"Input {missing_input_values} required for "
                        "plasma but not given when "
                        "instantiating the "
                        "plasma"
                    )
                current_property_object = plasma_property()
            else:
                current_property_object = plasma_property(self)
            for output in plasma_property.outputs:
                self.outputs_dict[output] = current_property_object
                plasma_property_objects.append(current_property_object)
        return plasma_property_objects

    def store_previous_properties(self):
        for property in self.previous_iteration_properties:
            base_property = property.outputs[0][9:]
            self.outputs_dict[property.outputs[0]].set_value(
                self.get_value(base_property)
            )

    def update(self, **kwargs):
        for key in kwargs:
            if key not in self.outputs_dict:
                raise PlasmaMissingModule(
                    f"Trying to update property {key} that is unavailable"
                )
            self.outputs_dict[key].set_value(kwargs[key])

        for module_name in self._resolve_update_list(kwargs.keys()):
            module = self.plasma_properties_dict[module_name].__class__
            if (
                not issubclass(module, ConvergedPlasmaProperty)
                or self.plasma_converged
            ):
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

        descendants_ob.sort(key=lambda val: sort_order.index(val))

        logger.debug(
            f"Updating modules in the following order:"
            f"{'->'.join(descendants_ob)}"
        )

        return descendants_ob

    def write_to_dot(self, fname, latex_label=True):
        #        self._update_module_type_str()

        try:
            import pygraphviz
        except:
            logger.warning(
                "pygraphviz missing. Plasma graph will not be generated."
            )
            return
        print_graph = self.graph.copy()
        print_graph = self.remove_hidden_properties(print_graph)
        for node in print_graph:
            print_graph.node[str(node)]["label"] = node
            if hasattr(self.plasma_properties_dict[node], "latex_formula"):
                formulae = self.plasma_properties_dict[node].latex_formula
                for output in range(len(formulae)):
                    formula = formulae[output]
                    label = formula.replace("\\", "\\\\")
                    print_graph.node[str(node)]["label"] += "\\n$"
                    print_graph.node[str(node)]["label"] += label
                    print_graph.node[str(node)]["label"] += "$"

        nx.write_dot(print_graph, fname)

    def write_to_tex(self, fname_graph):
        try:
            import dot2tex
        except:
            logger.warning(
                "dot2tex missing. Plasma graph will not be generated."
            )
            return

        temp_fname = tempfile.NamedTemporaryFile().name

        self.write_to_dot(temp_fname)

        dot_string = open(temp_fname).read()

        open(fname_graph, "w").write(dot2tex.dot2tex(dot_string, texmode="raw"))

        for line in fileinput.input(fname_graph, inplace=1):
            print(
                line.replace(
                    r"\documentclass{article}",
                    r"\documentclass[class=minimal,border=20pt]{standalone}",
                ),
            )

        for line in fileinput.input(fname_graph, inplace=1):
            print(
                line.replace(r"\enlargethispage{100cm}", ""),
            )

    def remove_hidden_properties(self, print_graph):
        for item in self.plasma_properties_dict.values():
            module = self.plasma_properties_dict[item.name].__class__
            if issubclass(module, HiddenPlasmaProperty):
                output = module.outputs[0]
                for value in self.plasma_properties_dict.keys():
                    if output in getattr(
                        self.plasma_properties_dict[value], "inputs", []
                    ):
                        for input in self.plasma_properties_dict[
                            item.name
                        ].inputs:
                            try:
                                position = self.outputs_dict[
                                    input
                                ].outputs.index(input)
                                label = self.outputs_dict[input].latex_name[
                                    position
                                ]
                                label = "$" + label + "$"
                                label = label.replace("\\", "\\\\")
                            except:
                                label = input.replace("_", "-")
                            self.graph.add_edge(
                                self.outputs_dict[input].name,
                                value,
                                label=label,
                            )
                print_graph.remove_node(str(item.name))
        return print_graph


class StandardPlasma(BasePlasma):
    def __init__(
        self,
        number_densities,
        atom_data,
        time_explosion,
        nlte_config=None,
        ionization_mode="lte",
        excitation_mode="lte",
        w=None,
        link_t_rad_t_electron=0.9,
    ):
        pass
