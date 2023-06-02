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
 
1) ``lte``, in which
:math:`J_\textrm{blue} = \textrm{Blackbody}(T_\textrm{rad})`
 
2) ``nebular`` in which
:math:`J_\textrm{blue} = W \times \textrm{Blackbody}(T_\textrm{rad})`
 
3) ``detailed`` in which the :math:`J_\textrm{blue}`
are calculated using an estimator (this is described in :doc:`../../../physics/montecarlo/estimators`).
 
TARDIS currently supports three different kinds of line interaction: ``scatter`` --- a resonance scattering implementation,
``macroatom`` --- the most complex form of line interaction described in :ref:`macroatom` and ``downbranch`` a simplified
version of ``macroatom`` in which only downward transitions are allowed (see :ref:`lineinteraction`).
 
Finally, ``w_epsilon`` describes the dilution factor to use to calculate :math:`J_\textrm{blue}` that are 0, which
causes problems with the code (so :math:`J_\textrm{blue}` are set to a very small number).

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