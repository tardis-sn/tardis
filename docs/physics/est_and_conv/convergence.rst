.. _convergence:

***********
Convergence
***********

As explained in :doc:`estimators.ipynb`, after each iteration the values for radiative temperature and dilution factor are updated by calling the ``advance_state`` method on a ``Simulation`` object. The goal of this is to eventually have the radiative temperature and dilution factor converge to a single value so that the steady-state plasma state can be determined. To ensure that the simulation converges, TARDIS employs additional convergence strategies. Currently, only one convergence strategy is available: damped convergence. This will be described in the following sections.


T_rad and W
-----------

As discussed :doc:`here <estimators.ipynb>`, TARDIS uses estimators to calculate estimated radiative temperatures (:math:`T_\mathrm{rad}`) and dilution factors (:math:`W`) in each cell. While TARDIS can then update the plasma state using the estimated values, there is a good chance that these estimated values would “overshoot” the true value we want to converge to (for example, if the current value of the dilution factor in some cell the dilution factor is .4, and the true steady-state value TARDIS wants to find is .45, there is a good chance that the estimated value will be greater than .45). This could make the simulation take longer to converge or, at worst, make it so the simulation does not converge at all. To account for this, users can set (in the :ref:`convergence section <conv-config>` of the monte carlo configuration) a "damping constant" for both the radiative temperature (:math:`d_{T_\mathrm{rad}}`) and the dilution factor (:math:`d_W`). When ``advance_state`` is called, these quantities update as follows:

.. math::
    T_\mathrm{rad}_{updated} = T_\mathrm{rad}_{current} + d_{T_\mathrm{rad}}(T_\mathrm{rad}_{estimated}-T_\mathrm{rad}_{current})
    
.. math::
    W_{updated} = W_{current} + d_W(W_{estimated}-W_{current})

This means, for example, if the damping constant is .5, the updated value is halfway between the current value and the estimated value. If the damping constant is .7, the updated value is 70% of the way between the current value and the estimated value, and so on. If the damping constant is 1, then the updated value is exactly the estimated value, and if the damping constant is zero, the value stays the same throughout the simulation and is not updated.


T_inner
-------


Convergence Information
-----------------------

During the simulation, information about the how :math:`T_\mathrm{rad}`, :math:`W`, and :math:`T_\mathrm{inner}` are updated as well as a comparison of the total output luminosity and the requested luminosity are logged at the INFO level (see :doc:`../../io/optional/logging_configuration.ipynb`) to give users a better idea of how the convergence process is working.

In addition, TARDIS allows for the displaying of convergence plots, which allows users to visualize the convergence process for :math:`T_\mathrm{rad}`, :math:`W`, :math:`T_\mathrm{inner}`, and the total luminosity of the supernova being modeled. For more information, see :doc:`../../io/visualization/convergence_plot.ipynb`.


Convergence Criteria
--------------------

TARDIS also allows users to stop the simulation if the simulation reaches a certain level of convergence. To enable this, users must set ``stop_if_converged=True`` in the :ref:`convergence section <conv-config>` of the monte carlo configuration. Also in the configuration, the quantities ``hold_iterations``, ``threshold``, and ``fraction`` are be specified to determine convergence as follows:

For the simulation to be considered to have converged, for ``hold_iterations`` successive iterations, the estimated values of :math:`T_\mathrm{rad}`, :math:`W`, and :math:`T_\mathrm{inner}` may differ from the previous value by a fraction of at most ``threshold`` in at least ``fraction`` fraction of the shells (for :math:`T_\mathrm{inner}`, since there is only one value, the ``fraction`` part does not apply). For example, if ``hold_iterations=3``, ``threshold=0.05`` for all three quantities, and ``fraction=.8``, the simulation will be considered to have converged if for 3 successive iterations if :math:`T_\mathrm{rad}` and :math:`W` change by at most 5% in at least 80% of the shells, *and* :math:`T_\mathrm{inner}` changes by at most 5%. See the :ref:`convergence section <conv-config>` of the monte carlo configuration for default values of these quantities.

.. note::

    ``hold_iterations`` and ``fraction`` are universal quantities, i.e. they are each a single value that applies to :math:`T_\mathrm{rad}` and :math:`W`, and for ``hold_iterations`` also :math:`T_\mathrm{inner}`. ``threshold``, on the other hand, is supplied for each quantity seperately, so for instance you could require :math:`T_\mathrm{rad}` to change by less than 1%, :math:`W` to change by less than 3%, and :math:`T_\mathrm{inner}` to change by less than 5% for convergence to be reached.
    
.. note::

    To determine convergence, we compare the estimated value, **not** the updated value (which is related to the estimated value via the damping constant), with the previous value. If :math:`T_\mathrm{inner}` is locked (see the previous section), the estimated value will still be calculated so convergence can be checked as usual.


Custom Convergence
------------------

The custom convergence strategy option is not currently implimented in TARDIS.
