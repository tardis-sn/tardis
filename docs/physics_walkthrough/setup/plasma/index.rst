.. _plasma:

******
Plasma
******


The role of the plasma module is to determine the ionisation and excitation states of the elements of the 
supernova ejecta, given the basic structure, including the elemental abundances, densities and radiation temperature.
After the calculation of the plasma state, the :math:`\tau_{\textrm{sobolev}}` values can be calculated.

The TARDIS plasma structure inherits from the `BasePlasma` class. The Type IIP workflow uses the equilibrium-backed
standard plasma assembled by `IIPPlasmaSolverFactory`.
for generating a plasma from the information provided by `model`. A variety of different plasmas can be generated
depending on the options selected in the plasma section of the TARDIS config. file. The options currently considered
by the Legacy Plasma when creating the plasma calculation structure include:

plasma:
 * ionization: lte/nebular
 * excitation: lte/dilute-lte
 * line_interaction_type: scatter/downbranch/macroatom
 * helium_treatment: dilute-lte/recomb-nlte
 * nlte: [can provide list of ion species to be treated in NLTE, as well as specifying the use of the coronal_approximation/classical_nebular settings.

The standard property graph uses these options to construct a map of the necessary plasma parameters and their dependencies
(using `NetworkX <https://networkx.github.io/>`_). For Type IIP simulations, hydrogen continuum populations, rate
coefficients, heating and cooling, and transport continuum state are calculated through the equilibrium solver components.
Before Monte Carlo continuum estimators are available, Type IIP continuum initialization uses dilute-LTE excitation and
nebular ionization with charge neutrality. Estimator-backed continuum updates use the Monte Carlo bound-free estimators
from the previous transport step.

The continuum population solve treats each ejecta shell independently. At a fixed electron density and Sobolev escape
probability, the level and ion rate matrices produce an element-normalized solution for every element, including the
fully stripped charge state. Ionization edges without configured continuum rates use the nebular balance
:math:`n_e N_{j+1} - \phi_j N_j = 0`. For estimator-backed updates, the level populations and Sobolev escape
probabilities are coupled to a bounded scalar charge-neutrality solve in each shell. The resulting level populations
update the escape probabilities until both population and escape-probability changes converge. Previous level populations
and escape probabilities initialize this iteration; previous ion and electron populations are not held fixed in the
converged equations. A final charge solve at the converged escape probabilities ensures that the public population and
Sobolev quantities describe the same state, following the iterative statistical-
equilibrium strategy of :cite:`Lucy2001` and the repeated transition-rate
back-substitution used in Monte Carlo line-transfer calculations
(:cite:`Lucy2003`).
Exact-zero nebular ionization factors follow the IIP input convention and are replaced by :math:`10^{-10}` times the
smallest positive factor before matrix assembly; this is not a fallback for an otherwise disconnected rate system.

Coupled continuum population calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For an element with atomic number :math:`Z` and number density :math:`N_Z`,
the coupled solver uses the element-normalized state

.. math::

    \boldsymbol{y}_Z = (p_{0,0}, \ldots, p_{0,L-1}, Y_1, \ldots, Y_Z),

where the :math:`p_{0,i}` are neutral-atom level fractions and :math:`Y_j` is
the fraction in charge state :math:`j`. The final state includes the fully
stripped :math:`Y_Z` state, even though it has no bound levels. The level
fractions in each level-bearing ion are summed to recover its total ion
population. Element conservation is imposed by

.. math::

    \sum_i p_{0,i} + \sum_{j=1}^{Z} Y_j = 1.

At fixed electron density and fixed Sobolev escape probabilities, every
element is solved with one linear rate matrix. Bound-bound rates combine
radiative rates scaled by the escape probability with electron-density-scaled
collisional rates. Bound-free edges use the photoionization and recombination
coefficients, including the Lucy detailed-balance factor exactly once. Edges
without configured continuum rates use the nebular balance
:math:`n_e Y_{j+1} - \phi_j Y_j = 0`.

The independent element solutions are coordinated through scalar charge
closure in each shell:

.. math::

    Q(n_e) = \sum_Z N_Z \sum_j jY_{Z,j}(n_e) - n_e = 0.

The root is bracketed between zero and the fully ionized electron density.
After charge closure, the solver recalculates the stimulated-emission factor,
Sobolev optical depth, and escape probability from the new level populations.
This back-substitution is repeated until both the population and escape-
probability updates converge, with deterministic damping if an update grows.
The final population solve is performed at the final escape probabilities, so
the public level, ion, electron-density, and Sobolev outputs describe one
coherent state. This scalar/Picard decomposition avoids the full linearization
approach of :cite:`AuerMihalas1969`; accelerated lambda iteration is likewise
not required because the Monte Carlo radiation estimators are held fixed during
the local plasma update (see :cite:`RybickiHummer1991`).

Initialization and diagnostics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When bound-free estimators are unavailable, initialization uses dilute-LTE
excitation together with nebular ionization and charge neutrality. Once a
complete set of Monte Carlo bound-free estimators is available, the current
estimators and current blue-line intensities determine the coupled update.
Previous populations and escape probabilities are warm starts only; they are
not mixed into the converged equations. Partial estimator or previous-state
inputs are errors rather than partial updates.

The hydrogen ground-state/Lyman-continuum suppression is an explicit policy of
the bound-free rate assembly and is independent of the iteration count. Exact
zero nebular ionization factors are regularized using the IIP convention above.
The solver rejects nonfinite or materially negative populations, singular or
disconnected rate systems, unbracketed charge residuals, and failed final
Sobolev back-substitution. Population residuals, element conservation, charge
closure, and the public Sobolev identities are checked before the state is
published.

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
If the necessary Python modules (`PyGraphviz <https://pygraphviz.github.io/>`_ and `dot2tex <https://dot2tex.readthedocs.io/en/latest/>`_) are available, TARDIS can output a .dot and a .tex file at the beginning of each run that can be compiled to produce a PDF image of the plasma module graph via the :code:`write_to_dot()` and :code:`write_to_tex()` functions, respectively. The nodes on this graph are the names of plasma properties, e.g. `Levels`, `TauSobolev`, `IonNumberDensity`, along with a list of outputs from those properties. These nodes are connected by edges linking them with the sources of their inputs. The .tex file contains the name of the input/output linking the properties (e.g. `levels`, :math:`\tau_{\textrm{sobolev}}`, :math:`n_{e}`), as well as the equations used to calculate them, written in LaTeX. See the tutorial below!

.. toctree::
    :maxdepth: 2

    ../../../how_to/output/how_to_plasma_graph

Updating the Plasma
-------------------
During each iteration of the main code, TARDIS updates the plasma using the `update_radiationfield` function. This requires, at minimum, new values for `t_rad` (the radiation temperature), `w` (the dilution factor) and `j_blues` (the intensity in the blue part of each line).

.. _plasma_calculations:

Plasma Calculations
-------------------

.. note::
    In this documentation we use the indices :math:`i, j, k` to mean atomic number, ion number and level number respectively.

`BasePlasma` serves as the base class for all plasmas and can just calculate the atom number densities for a given input of abundance fraction.

.. math::
    N_{atom} = \rho_\textrm{total} \times \textrm{Abundance fraction} / m_\textrm{atom}

In the next step the line and level tables are purged of entries that are not represented in the abundance fractions are saved in `BasePlasma.levels` and `BasePlasma.lines`. Finally, the function `BasePlasma.update_t_rad` is called at the end of initialization to update the plasma conditions to a new :math:`T_\textrm{radiation field}` (with the give t_rad). This function is the same in the other plasma classes and does the main part of the calculation. In the case of `BasePlasma` this is only setting `BasePlasma.beta_rad` to :math:`\frac{1}{k_\textrm{B}T_\textrm{rad}}`.

The next more complex class is `LTEPlasma` which will calculate the ionization balance and level populations in Local Thermal Equilibrium conditions (LTE). The :class:`NebularPlasma`-class inherits from `LTEPlasma` and uses a more complex description of the BasePlasma.

.. toctree::
    :maxdepth: 2

    lte_plasma
    nebular_plasma

TARDIS also allows for NLTE treatments of specified species, as well as special NLTE treatments for Helium.

.. note::
    The NLTE treatment of specified species is currently incompatible with the NLTE treatment for helium and cannot be used simultaneously.

.. toctree::
    :maxdepth: 2

    nlte
    helium_nlte
    two_photon_archive
    adiabatic_cooling_archive
    


.. _tau_sobolev:

Sobolev optical depth
---------------------

After the above calculations, TARDIS calculates the Sobolev optical depth :math:`\tau_\textrm{Sobolev}` with the following formula:


.. math::
    C_\textrm{Sobolev} = \frac{\pi e^2}{m_e c}

    \tau_\textrm{Sobolev} = C_\textrm{Sobolev}\,  \lambda\, f_{\textrm{lower}\rightarrow\textrm{upper}}\,
        t_\textrm{explosion}\, N_\textrm{lower}
        (1 - \frac{g_\textrm{lower}}{g_\textrm{upper}}\frac{N_\textrm{upper}}{N_\textrm{lower}})
        

Macro Atom Line Interaction Treatment
-------------------------------------

The following page describes the macro atom treatment of line interactions:

.. toctree::
    :maxdepth: 2
    
    macroatom
