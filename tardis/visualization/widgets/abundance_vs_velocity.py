import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go


def plot_abundance_vs_velocity(sim):
    """
    Plots the abundance of elements vs velocity provided simulation result using matplotlib.

    Parameters
    ----------
    sim : Model Object containing

        config : str or dict
            filename of configuration yaml file or dictionary

        atom_data : str or tardis.atomic.AtomData
            if atom_data is a string it is interpreted as a path to a file storing
            the atomic data. Atomic data to use for this TARDIS simulation.
            If set to None, the atomic data will be loaded according to keywords set in the configuration
            [default=None]


    Returns
    -------
    None

    """
    # code goes here

    # loading abundance sata in 'df'
    df = (sim.model.abundance).T

    # saving the velocity array in 'velocity'
    velocity = sim.model.velocity

    # obtaing a map of atomic numbers to element symbol
    # the index of a will be the atomic number and corresponding value will be element symbol
    a = atomic_number_to_element_symbol()

    # plotting the graph of abundance of element vs velocity

    # for loop iterate and plot abundance of each element at a time
    for index, row in df.T.iterrows():
        plt.plot(velocity[:-1], df[index], "o-", label=a[index])
    plt.legend()
    plt.title("Abundance of elements v/s Velocity")
    plt.xlabel("Velocity in " + str(sim.model.velocity.unit))
    plt.ylabel("Abundance (in Fraction)")
    plt.show()


def plotly_abundance_vs_velocity(sim):
    """
    Plots the abundance of elements vs velocity provided simulation result using plotly.

    Parameters
    ----------
    sim : Model Object containing

        config : str or dict
            filename of configuration yaml file or dictionary

        atom_data : str or tardis.atomic.AtomData
            if atom_data is a string it is interpreted as a path to a file storing
            the atomic data. Atomic data to use for this TARDIS simulation.
            If set to None, the atomic data will be loaded according to keywords set in the configuration
            [default=None]


    Returns
    -------
    None

    """

    # loading abundance sata in 'df'
    df = (sim.model.abundance).T

    # saving the velocity array in 'velocity'
    velocity = sim.model.velocity

    # obtaing a map of atomic numbers to element symbol
    # the index of a will be the atomic number and corresponding value will be element symbol
    a = atomic_number_to_element_symbol()

    # plotting the graph of abundance of element vs velocity

    # for loop iterate and plot abundance of each element at a time
    fig = go.Figure()
    for index, row in df.T.iterrows():
        fig.add_trace(
            go.Scatter(
                x=velocity[:-1],
                y=df[index],
                mode="lines+markers",
                name=a[index],
            )
        )
    fig.update_layout(
        xaxis=dict(
            title="Velocity in " + str(sim.model.velocity.unit),
            exponentformat="e",
        ),
        yaxis=dict(title="Abundance (in Fraction)"),
        title={
            "text": "Abundance of elements v/s Velocity",
            "y": 0.9,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        },
    )
    fig.update_yaxes(range=[0, 1])
    fig.show()


def atomic_number_to_element_symbol():
    """
    Convert atomic numbers into atomic symbols.

    Parameters
    ----------
    None


    Returns
    -------
    a : list
    the list whose index is mapped to atomic symbols.

    """
    a = [None] * 119
    a[1:16] = [
        "H",
        "He",
        "Li",
        "Be",
        "B",
        "C",
        "N",
        "O",
        "F",
        "Ne",
        "Na",
        "Mg",
        "Al",
        "Si",
        "P",
    ]
    a[16:31] = [
        "S",
        "Cl",
        "Ar",
        "K",
        "Ca",
        "Sc",
        "Ti",
        "V",
        "Cr",
        "Mn",
        "Fe",
        "Co",
        "Ni",
        "Cu",
        "Zn",
    ]
    a[31:45] = [
        "Ga",
        "Ge",
        "As",
        "Se",
        "Br",
        "Kr",
        "Rb",
        "Sr",
        "Y",
        "Zr",
        "Nb",
        "Mo",
        "Tc",
        "Ru",
    ]
    a[45:59] = [
        "Rh",
        "Pd",
        "Ag",
        "Cd",
        "In",
        "Sn",
        "Sb",
        "Te",
        "I",
        "Xe",
        "Cs",
        "Ba",
        "La",
        "Ce",
    ]
    a[59:73] = [
        "Pr",
        "Nd",
        "Pm",
        "Sm",
        "Eu",
        "Gd",
        "Tb",
        "Dy",
        "Ho",
        "Er",
        "Tm",
        "Yb",
        "Lu",
        "Hf",
    ]
    a[73:87] = [
        "Ta",
        "W",
        "Re",
        "Os",
        "Ir",
        "Pt",
        "Au",
        "Hg",
        "Tl",
        "Pb",
        "Bi",
        "Po",
        "At",
        "Rn",
    ]
    a[87:101] = [
        "Fr",
        "Ra",
        "Ac",
        "Th",
        "Pa",
        "U",
        "Np",
        "Pu",
        "Am",
        "Cm",
        "Bk",
        "Cf",
        "Es",
        "Fm",
    ]
    a[101:114] = [
        "Md",
        "No",
        "Lr",
        "Rf",
        "Db",
        "Sg",
        "Bh",
        "Hs",
        "Mt",
        "Ds",
        "Rg",
        "Cn",
        "Uut",
    ]
    a[114:119] = ["Fl", "Uup", "Lv", "Uus", "Uuo"]
    return a
