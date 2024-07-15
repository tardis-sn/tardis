import os
import re
import logging
import tempfile
import fileinput

import networkx as nx

from tardis.plasma.exceptions import PlasmaMissingModule, NotInitializedModule
from tardis.plasma.properties.base import *
from tardis.io.util import PlasmaWriterMixin

logger = logging.getLogger(__name__)


class BasePlasma(PlasmaWriterMixin):

    outputs_dict = {}
    hdf_name = "plasma"

    def __init__(self, plasma_properties, property_kwargs=None, **kwargs):
        self.outputs_dict = {}
        self.input_properties = []
        self.plasma_properties = self._init_properties(
            plasma_properties, property_kwargs, **kwargs
        )
        self._build_graph()
        self.update(**kwargs)

    def __getattr__(self, item):
        if item in self.outputs_dict:
            return self.get_value(item)
        else:
            super(BasePlasma, self).__getattribute__(item)

    def __setattr__(self, key, value):
        if key != "module_dict" and key in self.outputs_dict:
            raise AttributeError(
                "Plasma inputs can only be updated using " "the 'update' method"
            )
        else:
            super(BasePlasma, self).__setattr__(key, value)

    def __dir__(self):
        attrs = [item for item in self.__dict__ if not item.startswith("_")]
        attrs += [
            item for item in self.__class__.__dict__ if not item.startswith("_")
        ]
        attrs += self.outputs_dict.keys()
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
        # Adding all nodes
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
                        f" to this plasma"
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

    def _init_properties(
        self, plasma_properties, property_kwargs=None, **kwargs
    ):
        """
        Builds a dictionary with the plasma module names as keys

        Parameters
        ----------
        plasma_modules : list
            list of Plasma properties
        property_kwargs : dict
            dict of plasma module : kwargs pairs. kwargs should be a dict
            of arguments that will be passed to the __init__ method of
            the respective plasma module.
        kwargs : dictionary
            input values for input properties. For example,
            t_rad=[5000, 6000,], j_blues=[..]
        """
        if property_kwargs is None:
            property_kwargs = {}
        plasma_property_objects = []
        self.previous_iteration_properties = []
        self.outputs_dict = {}
        for plasma_property in plasma_properties:

            if issubclass(plasma_property, PreviousIterationProperty):
                current_property_object = plasma_property(
                    **property_kwargs.get(plasma_property, {})
                )
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
                        f"plasma but not given when "
                        f"instantiating the "
                        f"plasma"
                    )
                current_property_object = plasma_property(
                    **property_kwargs.get(plasma_property, {})
                )
            else:
                current_property_object = plasma_property(
                    self, **property_kwargs.get(plasma_property, {})
                )
            for output in plasma_property.outputs:
                self.outputs_dict[output] = current_property_object
                plasma_property_objects.append(current_property_object)
        return plasma_property_objects

    def store_previous_properties(self):
        for property in self.previous_iteration_properties:
            p = property.outputs[0]
            self.outputs_dict[p].set_value(
                self.get_value(re.sub(r"^previous_", "", p))
            )

    def update(self, **kwargs):
        for key in kwargs:
            if key not in self.outputs_dict:
                raise PlasmaMissingModule(
                    f"Trying to update property {key}" f" that is unavailable"
                )
            self.outputs_dict[key].set_value(kwargs[key])

        for module_name in self._resolve_update_list(kwargs.keys()):
            self.plasma_properties_dict[module_name].update()

    def freeze(self, *args):
        """
        Freeze plama properties.

        This method freezes plasma properties to prevent them from being
        updated: the values of a frozen property are fixed in the plasma
        calculation. This is useful for example for setting up test cases.

        Parameters
        ----------
        args : iterable of str
            Names of plasma properties to freeze.

        Examples
        --------
        >>> plasma.freeze('t_electrons')
        """
        for key in args:
            if key not in self.outputs_dict:
                raise PlasmaMissingModule(
                    "Trying to freeze property {0}"
                    " that is unavailable".format(key)
                )
            self.outputs_dict[key].frozen = True

    def thaw(self, *args):
        """
        Thaw plama properties.

        This method thaws (unfreezes) plasma properties allowing them to be
        updated again.

        Parameters
        ----------
        args : iterable of str
            Names of plasma properties to unfreeze.

        Examples
        --------
        >>> plasma.thaw('t_electrons')
        """
        for key in args:
            if key not in self.outputs_dict:
                raise PlasmaMissingModule(
                    "Trying to thaw property {0}"
                    " that is unavailable".format(key)
                )
            self.outputs_dict[key].frozen = False

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
        changed_modules : list
            all modules changed in the plasma

        Returns
        -------
            : list
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
            f'{"->".join(descendants_ob)}'
        )

        return descendants_ob

    def write_to_dot(self, fname, args=None, latex_label=True):
        """
        This method takes the NetworkX Graph generated from the _build_graph
        method, converts it into a DOT code, and saves it to a file

        Parameters
        ----------
        fname: str
            the name of the file the graph will be saved to
        args: list
            a list of optional settings for displaying the
            graph written in DOT format
        latex_label: boolean
            enables/disables writing LaTeX equations and
            edge labels into the file.
        """

        try:
            import pygraphviz
        except:
            logger.warn(
                "pygraphviz missing. Plasma graph will not be " "generated."
            )
            return
        print_graph = self.graph.copy()
        print_graph = self.remove_hidden_properties(print_graph)

        for node in print_graph:
            if latex_label == True:
                if hasattr(self.plasma_properties_dict[node], "latex_formula"):
                    print_graph.nodes[str(node)][
                        "label"
                    ] = f"\\\\textrm{{{node}: }}"
                    node_list = self.plasma_properties_dict[node]
                    formulae = node_list.latex_formula
                    for output in range(0, len(formulae)):
                        formula = formulae[output]
                        label = formula.replace("\\", "\\\\")
                        print_graph.nodes[str(node)]["label"] += label
                else:
                    print_graph.nodes[str(node)][
                        "label"
                    ] = f"\\\\textrm{{{node}}}"
            else:
                print_graph.nodes[str(node)]["label"] = node

        for edge in print_graph.edges:
            label = print_graph.edges[edge]["label"]
            print_graph.edges[edge]["label"] = " "
            if latex_label == True:
                print_graph.edges[edge]["texlbl"] = label

        nx.drawing.nx_agraph.write_dot(print_graph, fname)

        for line in fileinput.FileInput(fname, inplace=1):
            if latex_label == True:
                print(
                    line.replace(
                        r'node [label="\N"]',
                        r'node [texmode="math"]',
                    ),
                    end="",
                )
            else:
                print(
                    line.replace(
                        r'node [label="\N"];',
                        "",
                    ),
                    end="",
                )

        if args is not None:
            with open(fname, "r") as file:
                lines = file.readlines()

            for newline in args:
                lines.insert(1, f"\t{newline};\n")

            with open(fname, "w") as f:
                lines = "".join(lines)
                f.write(lines)

    def write_to_tex(self, fname_graph, scale=0.5, args=None, latex_label=True):
        """
        This method takes the NetworkX Graph generated from the _build_graph
        method, converts it into a LaTeX friendly format,
        and saves it to a file

        Parameters
        ----------
        fname_graph: str
            the name of the file the graph will be saved to
        args: list
            a list of optional settings for displaying the
            graph written in DOT format
        scale: float
            a scaling factor to expand/contract the generated
            graph
        latex_label: boolean
            enables/disables writing LaTeX equations and
            edge labels into the file.
        """
        try:
            import dot2tex
        except:
            logger.warn(
                "dot2tex missing. Plasma graph will not be " "generated."
            )
            return

        temp_fname = tempfile.NamedTemporaryFile().name

        self.write_to_dot(temp_fname, args=args, latex_label=latex_label)

        with open(temp_fname, "r") as file:
            dot_string = file.read().replace("\\\\", "\\")

        texcode = dot2tex.dot2tex(
            dot_string, format="tikz", crop=True, valignmode="dot"
        )

        with open(fname_graph, "w") as file:
            file.write(texcode)

        for line in fileinput.input(fname_graph, inplace=1):
            print(
                line.replace(
                    r"\documentclass{article}",
                    r"\documentclass[class=minimal,border=20pt]{standalone}",
                ),
                end="",
            )

        for line in fileinput.input(fname_graph, inplace=1):
            print(line.replace(r"\enlargethispage{100cm}", ""), end="")

        for line in fileinput.input(fname_graph, inplace=1):
            print(
                line.replace(
                    r"\begin{tikzpicture}[>=latex',line join=bevel,]",
                    r"\begin{tikzpicture}"
                    r"[>=latex',line join=bevel,"
                    rf"scale={scale}]",
                ),
                end="",
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
