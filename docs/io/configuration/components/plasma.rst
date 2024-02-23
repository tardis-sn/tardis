.. _plasma-config:

********************
Plasma Configuration
********************

The plasma configuration gives TARDIS the necessary information to calculate the plasma state (see :ref:`plasma`):

.. jsonschema:: schemas/plasma.yml

``inital_t_inner`` is initial temperature (will be updated for most modes of TARDIS --- see convergence section) of the black-body on the inner
boundary. ``initial_t_rad`` is the initial radiation temperature (will be updated for most modes of TARDIS - see convergence section). For debugging purposes and to compare to
:term:`synapps` calculations one can disable the electron scattering. TARDIS will issue a warning that this is not physical.
There are currently two ``plasma_type`` options available: ``nebular`` and ``lte``, which tell TARDIS how to run the
ionization equilibrium and level population calculations (see :ref:`plasma` for more information).
The radiative rates describe how to calculate the :math:`J_\textrm{blue}` needed for the :ref:`nlte` calculations and
:ref:`macroatom` calculations. There are three options for ``radiative_rates_type``: 
 
1) ``blackbody``, in which
:math:`J_\textrm{blue} = \textrm{Blackbody}(T_\textrm{rad})`
 
2) ``dilute-blackbody`` in which
:math:`J_\textrm{blue} = W \times \textrm{Blackbody}(T_\textrm{rad})`
 
3) ``detailed`` in which the :math:`J_\textrm{blue}`
are calculated using an estimator (this is described in :doc:`../../../physics/montecarlo/estimators`).
 
TARDIS currently supports three different kinds of line interaction: ``scatter`` --- a resonance scattering implementation,
``macroatom`` --- the most complex form of line interaction described in :ref:`macroatom` and ``downbranch`` a simplified
version of ``macroatom`` in which only downward transitions are allowed (see :ref:`lineinteraction`).
 
Finally, ``w_epsilon`` describes the dilution factor to use to calculate :math:`J_\textrm{blue}` that are 0, which
causes problems with the code (so :math:`J_\textrm{blue}` are set to a very small number).

Continuum Interaction
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

    plasma:
        link_t_rad_t_electron: 1.0
        continuum_interaction:
            species:
                - H I
                - H II
                - He I
                - He II 
            enable_adiabatic_cooling: True

This will add continuum interactions for all specified species. Setting :math:`T_\textrm{rad} = T_\textrm{electron}` through 
``link_t_rad_t_electron: 1.0`` is recommended to enforce LTE (unless the simulation uses NLTE treatment). 
``enable_adiabatic_cooling`` enables adiabatic cooling.

NLTE
^^^^

.. code-block:: yaml

    nlte:
        coronal_approximation: True
        classical_nebular: False

The NLTE configuration currently allows setting ``coronal_approximation``, which sets all :math:`J_\textrm{blue}` to 0.
This is useful for debugging with :term:`chianti` for example. Furthermore, one can enable 'classical_nebular' to set all
:math:`\beta_\textrm{Sobolev}` to 1. Both options are used for checking with other codes and should not be enabled in
normal operations.

NLTE Ionization
^^^^^^^^^^^^^^^

.. code-block:: yaml

    plasma:
        nlte_ionization_species: [H I, H II, He I, He II]
        nlte_solver: root
    
This option allows the user to specify which species should be included in the NLTE ionization treatment. Note that the
species must be present in the continuum interaction species as well.
Here, ``nlte_solver`` can be set to ``root`` or ``lu``. ``root`` is the default and uses a root solver to calculate the
NLTE populations. ``lu`` uses an iterative LU decomposition scheme to calculate the NLTE populations.

.. note ::

   ``lu`` iterates over the solutions up to a set tolerance. This tolerance is currently hard-coded to 1e-3. This
   can be changed in the code by changing the ``NLTE_POPULATION_SOLVER_TOLERANCE`` constant in ``tardis/plasma/properties/nlte_rate_equation_solver.py``.
   Furthermore, the maximum number of iterations is set to 1000. This can be changed in the code by changing the ``NLTE_POPULATION_SOLVER_MAX_ITERATIONS``
   constant in ``tardis/plasma/properties/nlte_rate_equation_solver.py``.

.. warning ::

    ``lu`` is generally faster than ``root`` but does not solve explicitly for the electron density. Therefore, it is
    not recommended to use ``lu`` for simulations where the electron density is important (e.g. for simulations where
    NLTE excitation is important).
