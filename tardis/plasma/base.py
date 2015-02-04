import networkx as nx

class BasePlasma(object):
        def __init__(self, plasma_properties):
            self._build_graph(plasma_properties)

        def _build_graph(self, plasma_properties):
            self.graph = nx.DiGraph()
            for plasma_property in plasma_properties:
                for input in plasma_property.inputs:
                    self.graph.add_edge(input, plasma_property.name)
            1/0


class PlasmaMissingModule(Exception):
    pass

class StandardPlasma(BasePlasma):

    def __init__(self, number_densities, atom_data, time_explosion,
                 delta_treatment=None, nlte_config=None, ionization_mode='lte',
                 excitation_mode='lte'):

        pass

