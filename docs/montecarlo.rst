Radiative Monte Carlo
=====================

.. :currentmodule:: tardis.montecarlo_multizone

The radiative monte carlo is initiated once the model is constructed.

Different line interactions


line_interaction_id == 0: scatter
line_interaction_id == 1: downbranch
line_interaction_id == 2: macro

Radiationfield estimators
-------------------------

During the monte-carlo run wie collect two estimators for the radiation field:

.. math::

    J_\textrm{estimator} &= \sum{\epsilon l}
    \bar{\nu}_\textrm{estimator} &=  \sum{\epsilon \nu l},

where :math:`\epsilon, \nu` are comoving energy and comoving frequency of a packet respectively.



These estimators help us t t