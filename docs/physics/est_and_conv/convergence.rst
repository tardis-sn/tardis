.. _convergence:

***********
Convergence
***********

:ref:`As mentioned <est_and_conv>`, one of the goals of TARDIS is to converge to a steady-state plasma state by using :doc:`monte carlo estimators <estimators.ipynb>`. However, to ensure that the simulation converges, TARDIS employs additional convergence strategies. Currently, only one convergence strategy is available: damped convergence. This will be described in the following sections. CONVERGENCE PLOTS????


T_rad and W
-----------

As discussed :doc:`here <estimators.ipynb>`, TARDIS uses estimators to calculate estimated radiative temperatures (:math:`T_\mathrm{rad}`) and dilution factors (:math:`W`) in each cell. While TARDIS can then update the plasma state using the estimated values, there is a good chance that these estimated values would “overshoot” the true value we want to converge to (for example, if the current value of the dilution factor in some cell the dilution factor is .4, and the true steady-state value TARDIS wants to find is .45, there is a good chance that the estimated value will be greater than .45). This could make the simulation take longer to converge or, at worst, make it so the simulation does not converge at all. To account for this, users can set (in the :ref:`convergence section of the monte carlo configuration <???>`) a damping constant???? for both the radiative temperature (:math:`d_{T_\mathrm{rad}}`) and the dilution factor (:math:`d_W`). When ``advance_state`` is called, these quantities update as follows:

.. math::
    T_\mathrm{rad}_{updated} = T_\mathrm{rad}_{current} + d_{T_\mathrm{rad}}(T_\mathrm{rad}_{estimated}-T_\mathrm{rad}_{current})
    
.. math::
    W_{updated} = W_{current} + d_W(W_{estimated}-W_{current})

This means, for example, if the damping constant is .5, the updated value is halfway between the current value and the estimated value. If the damping constant is .7, the updated value is 70% of the way between the current value and the estimated value, and so on. If the damping constant is 1, then the updated value is exactly the estimated value, and if the damping constant is zero, the value stays the same throughout the simulation and is not updated.


T_inner
-------

The temperature of the inner boundary is also updated throughout the simulation.

FINISH THIS!!! (mention lock cycles) (note about the exponent being negative)

!!!!Mention J_blues in a note?!!!!


Convergence Criteria
--------------------

TARDIS also allows users to stop the simulation if the simulation reaches a certain level of convergence. To enable this, users must set ________. In the convergence section ………. users can set a threshold for each of the three relevant quantities (tradtinnerw………). If the value of the dilution factor changes by a ESTIMATED VALUE?????

Users can also set what fraction of cells must the radiative temperature and dilution factor converge in for the simulation to be considered to have converged (i.e. if this is set to 0.8, then for the simulation to be considered to have converge, 

all of this is checked during advance state


Custom Convergence
------------------

The custom convergence strategy option is not currently implimented in TARDIS.



USE THIS???:
As explained in :doc:`estimators.ipynb`, after each iteration the values for radiative temperature and dilution factor are updated by calling the ``advance_state`` method on a ``Simulation`` object. The goal of this is to eventually have the radiative temperature and dilution factor converge to a single value.
