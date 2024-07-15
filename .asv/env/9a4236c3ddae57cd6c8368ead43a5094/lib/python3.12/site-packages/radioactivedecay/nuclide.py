"""
The nuclide module defines the ``Nuclide`` class. Each ``Nuclide`` instance
represents one nuclide, and can be contructed from a nuclide string or
`zzzaaassss` canonical id. The class properties provide an easy way of
extracting atomic number and mass information, as well as a clean name string.
Furthermore, additional methods provide an access point for mass data of
nuclides, and the decay data of radionuclides, if present in a specified
dataset. The default decay dataset used if none is supplied to the constructor
is rd.DEFAULTDATA.

The docstring code examples assume that ``radioactivedecay`` has been imported
as ``rd``:

.. highlight:: python
.. code-block:: python

    >>> import radioactivedecay as rd

"""

from collections import deque
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib
import networkx as nx
import numpy as np

from radioactivedecay.decaydata import DEFAULTDATA, DecayData
from radioactivedecay.plots import (
    _check_fig_axes,
    _parse_decay_mode_label,
    _parse_nuclide_label,
)
from radioactivedecay.utils import build_id, elem_to_Z, parse_nuclide


class Nuclide:
    """
    ``Nuclide`` instances serve as name and atomic number/mass parsing
    objects for any nuclide or element.

    Parameters
    ----------
    nuclide : str or int
        Specify the nuclide with either a name string (e.g. 'H-3', 'H3' or '3H') or a canonical id
        (zzzaaassss format).
    decay_data : DecayData, optional
        Decay dataset (default is the ICRP-107 dataset).

    Attributes
    ----------
    nuclide: str
        Nuclide name string.
    decay_data : DecayData
        Decay dataset.

    Examples
    --------
    >>> rd.Nuclide('K-40')
    Nuclide: K-40, decay dataset: icrp107_ame2020_nubase2020
    >>> rd.Nuclide('K40')
    Nuclide: K-40, decay dataset: icrp107_ame2020_nubase2020
    >>> rd.Nuclide(190400000)
    Nuclide: K-40, decay dataset: icrp107_ame2020_nubase2020
    >>> rd.Nuclide(280560001)
    Nuclide: Ni-56m, decay dataset: icrp107_ame2020_nubase2020

    """

    def __init__(
        self, nuclide: Union[str, int], decay_data: DecayData = DEFAULTDATA
    ) -> None:
        self.decay_data = decay_data
        self.nuclide = parse_nuclide(
            nuclide, self.decay_data.nuclides, self.decay_data.dataset_name
        )

    @property
    def Z(self) -> int:
        """
        Returns the atomic number of the nuclide.

        Returns
        -------
        int
            Atomic number of the nuclide.

        Examples
        --------
        >>> H3 = rd.Nuclide('H-3')
        >>> H3.Z
        1

        """

        return elem_to_Z(self.nuclide.split("-")[0])

    @property
    def A(self) -> int:
        """
        Returns the mass number of the nuclide.

        Returns
        -------
        int
            Mass number of the nuclide.

        Examples
        --------
        >>> H3 = rd.Nuclide('H-3')
        >>> H3.A
        3

        """

        return int(self.nuclide.split("-")[1].strip("mn"))

    @property
    def state(self) -> str:
        """
        Returns the metastable state character, i.e. '' for ground state, 'm' for first metastable
        state, 'n' for second metastable state, etc..

        Returns
        -------
        str
            Metastable state character.

        Examples
        --------
        >>> H3 = rd.Nuclide('H-3')
        >>> H3.state
        ''
        >>> Ba137m = rd.Nuclide('Ba-137m')
        >>> Ba137m.state
        'm'

        """

        return self.nuclide.split("-")[1].strip("0123456789")

    @property
    def id(self) -> int:
        """
        Returns the canonical nuclide id, in zzzaaassss form. Ground state is 0000, first excited
        state ("m") is 0001, second ("n") is 0002, etc.

        Returns
        -------
        int
            Canonical id of the nuclide.

        Examples
        --------
        >>> H3 = rd.Nuclide('H-3')
        >>> H3.id
        10030000
        >>> Ba137m = rd.Nuclide('Ba-137m')
        >>> Ba137m.id
        561370001

        """

        return build_id(self.Z, self.A, self.state)

    @property
    def atomic_mass(self) -> float:
        """
        Returns the atomic mass of the nuclide, in g/mol.

        Returns
        -------
        float
            Atomic mass of the nuclide in g/mol.

        Examples
        --------
        >>> H3 = rd.Nuclide('H-3')
        >>> H3.atomic_mass
        3.01604928132
        >>> Ba137m = rd.Nuclide('Ba-137m')
        >>> Ba137m.atomic_mass
        136.9065375271172

        """

        return self.decay_data.scipy_data.atomic_masses[
            self.decay_data.nuclide_dict[self.nuclide]
        ]

    def half_life(self, units: str = "s") -> Union[float, str]:
        """
        Returns the half-life of a nuclide as a float in your chosen
        units, or as a human-readable string with appropriate units.

        Parameters
        ----------
        units : str, optional
            Units for half-life. Options are 'ps', 'ns', 'μs', 'us',
            'ms', 's', 'm', 'h', 'd', 'y', 'ky', 'My', 'By', 'Gy',
            'Ty', 'Py', and common spelling variations. Default is 's',
            i.e. seconds. Use 'readable' to get a string of the
            half-life in human-readable units.

        Returns
        -------
        float or str
            Half-life of the nuclide.

        Examples
        --------
        >>> K40 = rd.Nuclide('K-40')
        >>> K40.half_life('y')
        1251000000.0
        >>> K40.half_life('readable')
        '1.251 By'
        >>> Fe56 = rd.Nuclide('Fe-56')
        >>> Fe56.half_life('readable')
        'stable'

        """

        return self.decay_data.half_life(self.nuclide, units)

    def progeny(self) -> List[str]:
        """
        Returns list of the direct progeny of the nuclide.

        Returns
        -------
        list
            List of the direct progeny of the nuclide, ordered by decreasing branching fraction.

        Examples
        --------
        >>> K40 = rd.Nuclide('K-40')
        >>> K40.progeny()
        ['Ca-40', 'Ar-40']

        """

        return self.decay_data.progeny[self.decay_data.nuclide_dict[self.nuclide]]

    def branching_fractions(self) -> List[float]:
        """
        Returns the branching fractions to the direct progeny of a
        radionuclide.

        Returns
        -------
        list
            List of branching fractions.

        Examples
        --------
        >>> K40 = rd.Nuclide('K-40')
        >>> K40.branching_fractions()
        [0.8914, 0.1086]

        """

        return self.decay_data.bfs[self.decay_data.nuclide_dict[self.nuclide]]

    def decay_modes(self) -> List[str]:
        """
        Returns the decay modes for a radionuclide, as defined in the
        decay dataset. Note: the decay mode strings returned are not
        lists of all the different radiation types emitted during the
        parent to progeny decay processes. They are the labels defined
        in the decay dataset to classify the parent to progeny decay
        type (e.g. '\u03b1', '\u03b2-' or 'IT').

        Returns
        -------
        list
            List of decay modes.

        Examples
        --------
        >>> K40 = rd.Nuclide('K-40')
        >>> K40.decay_modes()
        ['\u03b2-', '\u03b2+ & EC']

        """

        return self.decay_data.modes[self.decay_data.nuclide_dict[self.nuclide]]

    def plot(
        self,
        label_pos: float = 0.5,
        fig: Optional[matplotlib.figure.Figure] = None,
        axes: Optional[matplotlib.axes.Axes] = None,
        kwargs_draw: Optional[Dict[str, Any]] = None,
        kwargs_edge_labels: Optional[Dict[str, Any]] = None,
    ) -> Tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
        """
        Plots a diagram of the decay chain of a radionuclide. Then
        creates a NetworkX DiGraph and plot of it using NetworkX's
        Matplotlib-based plotting functionality.

        Some of the NetworkX default plotting parameters are changed to
        produce nice decay chain diagrams. However, users retain
        control over these parameters via kwargs_draw and
        kwargs_edge_labels. For more information on the various
        NetworkX plotting parameters, refer to its `documentation
        <https://networkx.org/documentation/stable/reference/drawing.html>`_.

        Parameters
        ----------
        label_pos : float, optional
            Position of labels along edges. Default is 0.5. If you find
            that edge labels are overlapping in the decay chain
            diagram, try increasing this parameter to e.g. 0.66.
        fig : None or matplotlib.figure.Figure, optional
            matplotlib figure object to use, or None makes
            ``radioactivedecay`` create one (default is None).
        axes : None or matplotlib.axes.Axes, optional
            matplotlib axes object to use, or None makes
            ``radioactivedecay`` create one (default is None).
        **kwargs_draw, optional
            Keyword arguments for networkx.draw().
        **kwargs_edge_labels, optional
            Keyword arguments for networkx.draw_networkx_edge_labels().

        Returns
        -------
        fig : matplotlib.figure.Figure
            matplotlib figure object used to plot the decay chain.
        ax : matplotlib.axes.Axes
            matplotlib axes object used to plot the decay chain.

        """

        digraph, max_generation, max_xpos = _build_decay_digraph(self, nx.DiGraph())

        positions = nx.get_node_attributes(digraph, "pos")
        node_labels = nx.get_node_attributes(digraph, "label")
        edge_labels = nx.get_edge_attributes(digraph, "label")

        fig, axes = _check_fig_axes(
            fig, axes, figsize=(3 * max_xpos + 1.5, 3 * max_generation + 1.5)
        )

        if kwargs_draw is None:
            kwargs_draw = {}
        if "node_size" not in kwargs_draw:
            kwargs_draw["node_size"] = 6000
        if "node_color" not in kwargs_draw:
            kwargs_draw["node_color"] = "#FFFFFF"
        if "edgecolors" not in kwargs_draw:
            kwargs_draw["edgecolors"] = "#000000"

        nx.draw(
            G=digraph,
            pos=positions,
            ax=axes,
            labels=node_labels,
            **kwargs_draw,
        )

        if kwargs_edge_labels is None:
            kwargs_edge_labels = {}
        if "font_size" not in kwargs_edge_labels:
            kwargs_edge_labels["font_size"] = 12
        if "bbox" not in kwargs_edge_labels:
            kwargs_edge_labels["bbox"] = {
                "boxstyle": None,
                "ec": (1.0, 1.0, 1.0),
                "fc": (1.0, 1.0, 1.0),
            }
        if "rotate" not in kwargs_edge_labels:
            kwargs_edge_labels["rotate"] = False

        nx.draw_networkx_edge_labels(
            G=digraph,
            pos=positions,
            edge_labels=edge_labels,
            label_pos=label_pos,
            ax=axes,
            **kwargs_edge_labels,
        )

        axes.set_xlim(-0.3, max_xpos + 0.3)
        axes.set_ylim(-max_generation - 0.3, 0.3)

        return fig, axes

    def __repr__(self) -> str:
        rep = f"Nuclide: {self.nuclide}, decay dataset: {self.decay_data.dataset_name}"

        return rep

    def __eq__(self, other: object) -> bool:
        """
        Check whether two ``Nuclide`` instances are equal with ``==`` operator.
        """

        if not isinstance(other, Nuclide):
            return NotImplemented
        return self.nuclide == other.nuclide and self.decay_data == other.decay_data

    def __ne__(self, other: object) -> bool:
        """
        Check whether two ``Nuclide`` instances are not equal with ``!=`` operator.
        """

        if not isinstance(other, Nuclide):
            return NotImplemented
        return not self.__eq__(other)

    def __hash__(self) -> int:
        """
        Hash function for ``Nuclide`` instances.
        """

        return hash((self.nuclide, self.decay_data.dataset_name))


def _build_decay_digraph(
    parent: Nuclide,
    digraph: nx.classes.digraph.DiGraph,
) -> nx.classes.digraph.DiGraph:
    """
    Build a networkx DiGraph for the decay chain of this nuclide.

    Parameters
    ----------
    parent : Nuclide
        Nuclide instance of the parent nuclide of the decay chain.
    digraph : networkx.classes.digraph.DiGraph
        DiGraph for the decay chain.

    Returns
    -------
    digraph : networkx.classes.digraph.DiGraph
        DiGraph of the decay chain.
    max_generation : int
        Number of generations of progeny in the decay chain.
    max_xpos : int
        Maximum number of progeny within any one generation of the decay chain.

    """

    generation_max_xpos = {0: 0}

    dequeue = deque([parent.nuclide])
    generations = deque([0])
    xpositions = deque([0])
    node_label = (
        f"{_parse_nuclide_label(parent.nuclide)}\n{parent.half_life('readable')}"
    )
    digraph.add_node(parent.nuclide, generation=0, xpos=0, label=node_label)
    seen = {parent.nuclide}

    while len(dequeue) > 0:
        parent_name = dequeue.popleft()
        generation = generations.popleft() + 1
        xpos = xpositions.popleft()
        if generation not in generation_max_xpos:
            generation_max_xpos[generation] = -1
        parent = Nuclide(parent_name, parent.decay_data)

        progeny = parent.progeny()
        branching_fractions = parent.branching_fractions()
        decay_modes = parent.decay_modes()

        xpos = max(xpos, generation_max_xpos[generation] + 1)
        xcounter = 0
        for idx, prog in enumerate(progeny):
            if prog not in seen:
                node_label = _parse_nuclide_label(prog)
                if prog in parent.decay_data.nuclide_dict:
                    node_label += f"\n{parent.decay_data.half_life(prog, 'readable')}"
                    if np.isfinite(parent.decay_data.half_life(prog)):
                        dequeue.append(prog)
                        generations.append(generation)
                        xpositions.append(xpos + xcounter)
                if prog == "SF":
                    prog = parent.nuclide + "_SF"

                digraph.add_node(
                    prog,
                    generation=generation,
                    xpos=xpos + xcounter,
                    label=node_label,
                )
                seen.add(prog)

                if xpos + xcounter > generation_max_xpos[generation]:
                    generation_max_xpos[generation] = xpos + xcounter
                xcounter += 1

            edge_label = (
                _parse_decay_mode_label(decay_modes[idx])
                + "\n"
                + str(branching_fractions[idx])
            )
            digraph.add_edge(parent.nuclide, prog, label=edge_label)

    for node in digraph:
        digraph.nodes[node]["pos"] = (
            digraph.nodes[node]["xpos"],
            digraph.nodes[node]["generation"] * -1,
        )

    return digraph, max(generation_max_xpos), max(generation_max_xpos.values())
