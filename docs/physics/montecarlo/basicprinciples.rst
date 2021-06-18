.. _montecarlo_basics:

*************************************************
Monte Carlo Radiative Transfer - Basic Principles
*************************************************

Radiative transfer describes how the properties of the electromagnetic
radiation field change as it propagates through an ambient material. Due to the
complex coupling induced by the many interactions the photons which constitute
the radiation field may perform with the surrounding material, solving the
radiative transfer problem is in most astrophysical cases only possible by
means of numerical calculations. While there are many different numerical
techniques available to tackle this, Monte Carlo techniques have become a
successful and elegant tool particularly for radiative transfer problems in
supernovae.

Monte Carlo Radiative Transfer methods are probabilistic techniques and draw
inspiration from the microscopical interpretation of how photons propagate
through and interact with the ambient material. A sufficiently large number of
representative "machine photons" are considered and their propagation history
solved in a stochastic process. The initial properties of these photons are
randomly (in a probabilistic sense) assigned in accordance with the macroscopic
properties of the radiation field (see :ref:`Energy Packets <initialization>`)
and in a similar manner the decisions about when, where and how the machine
photons interact with the surrounding material are made (see :ref:`Propagation
<propagation>`). If this process is repeated for a large enough number of machine
photons, the ensemble behaviour and thus the macroscopic evolution of the
radiation field is recovered (see :ref:`Estimators <estimators>`).
