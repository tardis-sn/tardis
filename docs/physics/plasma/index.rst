.. _plasma:

******
Plasma
******


The role of the plasma module is to determine the ionisation and excitation states of the elements of the 
supernova ejecta, given the basic structure, including the elemental abundances, densities and radiation temperature.
After the calculation of the plasma state, the :math:`\tau_{\textrm{sobolev}}` values can be calculated.

The TARDIS plasma structure inherits from the `BasePlasma` class. The code currently uses the `LegacyPlasmaArray`
for generating a plasma from the information provided by `model`. A variety of different plasmas can be generated
depending on the options selected in the plasma section of the TARDIS config. file. The options currently considered
by the Legacy Plasma when creating the plasma calculation structure include:

plasma:
 * ionization: lte/nebular
 * excitation: lte/dilute-lte
 * line_interaction_type: scatter/downbranch/macroatom
 * helium_treatment: dilute-lte/recomb-nlte
 * nlte: [can provide list of ion species to be treated in NLTE, as well as specifying the use of the coronal_approximation/classical_nebular settings.

`LegacyPlasmaArray` uses these options to construct a map of the necessary plasma parameters that demonstrates how these parameters are dependent on one another (using `NetworkX <https://networkx.github.io/>`_). Each time a particular parameter of the plasma is updated, all of the parameters dependent (directly or indirectly) on that particular one can be easily updated automatically, without requiring that all the plasma calculations are repeated.

Properties, Inputs and Outputs
------------------------------
Each TARDIS plasma possesses an array of plasma properties, which are used to calculate plasma parameter values. Most plasma properties have a single output, e.g.
 * `GElectron`: (`g_electron`,)
 * `HeliumNLTE`: (`helium_population`,)

but some have two or more, e.g.
 * `IonNumberDensity`: (`ion_number_density`, `electron_densities`)
 * `Levels`: (`levels`, `excitation_energy`, `metastability`, `g`)

Every property has a `calculate` function that returns the values of its outputs. The arguments required for that function become the property inputs. TARDIS will raise an error if it does not have all of the required inputs for a particular property. It will also raise an error if there is an output loop, i.e. if two properties are dependent on each other. Some different properties share output names; for example, `PhiSahaLTE` and `PhiSahaNebular` both have an output called `phi`. That is because the `phi` value is calculated differently depending on the ionization method selected, but once calculated, both values interact in the same way with the rest of the plasma. TARDIS will import only one of the `phi` properties when initialising the plasma.

The Plasma Graph
----------------
If the necessary Python modules (`PyGraphviz <https://pygraphviz.github.io/>`_ and `dot2tex <https://dot2tex.readthedocs.io/en/latest/>`_) are available, TARDIS will automatically output a .tex file at the beginning of each run that can be compiled to produce a PDF image of the plasma module graph. The nodes on this graph are the names of plasma properties, e.g. `Levels`, `TauSobolev`, `IonNumberDensity`, along with a list of outputs from those properties and equations showing how they are calculated. These nodes are connected by arrows linking nodes with the sources of their inputs, and labelled with the name of the input/output linking the two properties, e.g. `levels`, :math:`\tau_{\textrm{sobolev}}`, :math:`n_{e}`.

Updating the Plasma
-------------------
During each iteration of the main code, TARDIS updates the plasma using the `update_radiationfield` function. This requires, at minimum, new values for `t_rad` (the radiation temperature), `w` (the dilution factor) and `j_blues` (the intensity in the blue part of each line).

.. _plasma_calculations:

Plasma Calculations
-------------------
.. toctree::
    :maxdepth: 2

    lte_plasma
    nebular_plasma
    macroatom
    nlte
    helium_nlte