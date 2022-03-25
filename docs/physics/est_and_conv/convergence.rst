.. _convergence:

***********
Convergence
***********

As explained in :doc:`estimators`, after each iteration the values for radiative temperature and dilution factor are updated by calling the ``advance_state`` method on a ``Simulation`` object. The goal of this is to eventually have the radiative temperature and dilution factor converge to a single value so that the steady-state plasma state can be determined. To ensure that the simulation converges, TARDIS employs additional convergence strategies. Currently, only one convergence strategy is available: damped convergence. This will be described in the following sections.

.. note::
    
    Unless otherwise noted, all user-supplied quantities mentioned on this page are supplied in the :ref:`convergence section of the Monte Carlo configuration <conv-config>`, which will be referenced as the convergence configuration.


T_rad and W
-----------

As discussed :doc:`here <estimators>`, TARDIS uses estimators to calculate estimated radiative temperatures (:math:`T_\mathrm{rad}`) and dilution factors (:math:`W`) in each cell. While TARDIS can then update the plasma state using the estimated values, there is a good chance that these estimated values would “overshoot” the true value we want to converge to (for example, if the current value of the dilution factor in some cell the dilution factor is .4, and the true steady-state value TARDIS wants to find is .45, there is a good chance that the estimated value will be greater than .45). This could make the simulation take longer to converge or, at worst, make it so the simulation does not converge at all. To account for this, users can set (in the convergence configuration) a "damping constant" for both the radiative temperature (:math:`d_{T_\mathrm{rad}}`) and the dilution factor (:math:`d_W`). When ``advance_state`` is called, these quantities update as follows:

.. math::
    T_\mathrm{rad\ updated} = T_\mathrm{rad\ current} + d_{T_\mathrm{rad}}(T_\mathrm{rad\ estimated}-T_\mathrm{rad\ current})
    
and
    
.. math::
    W_\mathrm{updated} = W_\mathrm{current} + d_W(W_\mathrm{estimated}-W_\mathrm{current}).

This means, for example, if the damping constant is .5, the updated value is halfway between the current value and the estimated value. If the damping constant is .7, the updated value is 70% of the way between the current value and the estimated value, and so on. **If the damping constant is 1, then the updated value is exactly the estimated value, and if the damping constant is zero, the value stays the same throughout the simulation and is not updated.**


T_inner
-------

The temperature of the inner boundary, :math:`T_\mathrm{inner}`, plays a unique role in the simulation, as it is the primary determiner of the output luminosity. This is because the the luminosity of the inner boundary is proportional to :math:`T_\mathrm{inner}^4` (see :doc:`../montecarlo/initialization`). Thus, :math:`T_\mathrm{inner}` is updated throughout the simulation in order to match the output luminosity to the requested luminosity specified in the :doc:`supernova configuration <../../io/configuration/components/supernova>` between the bounds specified in the supernova configuration. However, there is not necessarily a quartic relationship between :math:`T_\mathrm{inner}` and the output luminosity, as changing :math:`T_\mathrm{inner}` also changes the frequency distribution of the initialized packets (once again see :doc:`../montecarlo/initialization`). This then affects the light-matter interactions, affecting which packets make it to the outer boundary, which also affects the output luminosity. Because of this, there is not an exact way to estimate :math:`T_\mathrm{inner}`. To do this estimation, we use

.. math::
    T_\mathrm{inner\ estimated} = T_\mathrm{inner\ current} * \left(\frac{L_\mathrm{output}}{L_\mathrm{requested}}\right)^{\mathrm{t\_inner\_update\_exponent}}
    
where :math:`L_\mathrm{output}` is the output luminosity calculated by adding up the luminosity of each packet (see :doc:`../spectrum/basic`) between the bounds specified in the :doc:`supernova configuration <../../io/configuration/components/supernova>`, :math:`L_\mathrm{requested}` is the luminosity requested also in the supernova configuration (requested between those bounds previously mentioned), and ``t_inner_update_exponent`` is provided by the user in the convergence configuration. Note that what we are doing is "correcting" the previous value of the inner temperature by a factor of :math:`\left(\frac{L_\mathrm{output}}{L_\mathrm{requested}}\right)^{\mathrm{t\_inner\_update\_exponent}}`. Note that if :math:`\frac{L_\mathrm{output}}{L_\mathrm{requested}}` is greater than 1, we want to lower :math:`T_\mathrm{inner}` as the output luminosity is too high, and vice versa if the ratio is less than 1. Thus ``t_inner_update_exponent`` should be negative. Naively one might set ``t_inner_update_exponent=-0.25``, however as a default TARDIS uses ``t_inner_update_exponent=-0.5`` as -0.25 may undershoot the correct :math:`T_\mathrm{inner}` because of its previously mentioned effects on the initial frequency distribution.

After calculating the estimated :math:`T_\mathrm{inner}`, the quantity is updated using damped convergence with its own damping constant (once again set in the convergence configuration):

.. math::
    T_\mathrm{inner\ updated} = T_\mathrm{inner\ current} + d_{T_\mathrm{inner}}(T_\mathrm{inner\ estimated}-T_\mathrm{inner\ current}).

Once again, If the damping constant is 1, then the updated value is exactly the estimated value, and if the damping constant is zero, the value stays the same throughout the simulation and is not updated.

Additionally, because of the vast impact of :math:`T_\mathrm{inner}` on the simulation, one may want to update it less frequently -- i.e. allow :math:`W` and :math:`T_\mathrm{rad}` to reach a steady-state value for a particular :math:`T_\mathrm{inner}` before updating :math:`T_\mathrm{inner}`. To do this, in the convergence configuration we set ``lock_t_inner_cycles``, which is the number of iterations to wait before updating :math:`T_\mathrm{inner}`. It is set to 1 by default, meaning :math:`T_\mathrm{inner}` would be updated every iteration.


Convergence Information
-----------------------

During the simulation, information about the how :math:`T_\mathrm{rad}`, :math:`W`, and :math:`T_\mathrm{inner}` are updated as well as a comparison of the total output luminosity and the requested luminosity are logged at the INFO level (see :doc:`../../io/optional/logging_configuration`) to give users a better idea of how the convergence process is working.

In addition, TARDIS allows for the displaying of convergence plots, which allows users to visualize the convergence process for :math:`T_\mathrm{rad}`, :math:`W`, :math:`T_\mathrm{inner}`, and the total luminosity of the supernova being modeled. For more information, see :doc:`../../io/visualization/convergence_plot`.


Convergence Criteria
--------------------

TARDIS also allows users to stop the simulation if the simulation reaches a certain level of convergence. To enable this, users must set ``stop_if_converged=True`` in the convergence configuration. Also in the configuration, the quantities ``hold_iterations``, ``threshold``, and ``fraction`` are be specified to determine convergence as follows:

For the simulation to be considered to have converged, for ``hold_iterations`` successive iterations, the estimated values of :math:`T_\mathrm{rad}`, :math:`W`, and :math:`T_\mathrm{inner}` may differ from the previous value by a fraction of at most ``threshold`` in at least ``fraction`` fraction of the shells (for :math:`T_\mathrm{inner}`, since there is only one value, the ``fraction`` part does not apply). For example, if ``hold_iterations=3``, ``threshold=0.05`` for all three quantities, and ``fraction=.8``, the simulation will be considered to have converged if for 3 successive iterations the estimated values of :math:`T_\mathrm{rad}` and :math:`W` differ from the current respective values by at most 5% in at least 80% of the shells, *and* the estimated :math:`T_\mathrm{inner}` differs by at most 5%. See the :ref:`convergence configuration schema <conv-config>` for default values of these quantities.

.. note::

    To determine convergence, we compare the estimated value, **not** the updated value (which is related to the estimated value via the damping constant), with the previous value. If :math:`T_\mathrm{inner}` is locked (see the previous section), the estimated value will still be calculated so convergence can be checked as usual.


.. note::

    ``hold_iterations`` and ``fraction`` are universal quantities, i.e. they are each a single value that applies to :math:`T_\mathrm{rad}` and :math:`W`, and for ``hold_iterations`` also :math:`T_\mathrm{inner}`. ``threshold``, on the other hand, is supplied for each quantity separately, so for instance you could require :math:`T_\mathrm{rad}` to differ by less than 1%, :math:`W` to differ by less than 3%, and :math:`T_\mathrm{inner}` to differ by less than 5% for convergence to be reached.
    

Custom Convergence
------------------

The custom convergence strategy option is not currently implemented in TARDIS.
