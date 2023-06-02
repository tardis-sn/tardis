.. _montecarlo_basics:

*************************************************
Monte Carlo Radiative Transfer - Basic Principles
*************************************************

Radiative transfer describes how electromagnetic radiation (light) propagates through a medium. Since there are
a large number of light-matter interactions that have the possibility of affecting the light propagation, solving
the radiative transfer problem is in most astrophysical cases only possible by means of numerical calculations.
While there are many different numerical techniques available to tackle this, Monte Carlo techniques have become a
successful and elegant tool particularly for radiative transfer problems in supernovae.

Monte Carlo Radiative Transfer methods track a sufficiently large number of photons (light particles) as they
propagate through the supernova ejecta. The initial properties of these photons are randomly (in a probabilistic
sense) assigned in accordance with the macroscopic properties of the radiation field (see :doc:`initialization`)
and in a similar manner the decisions about when, where and how the photons interact with the surrounding material
are made (see :ref:`Propagation <propagation>`). Given a large enough sample, these photons behave as a microcosom
of all of the light traveling through the ejecta -- that is, based on the behavior of these photons, we can draw
conclusions about the propagation of light through the ejecta as a whole (see :ref:`estimators`). This is eventually
used to determine the actual steady-state plasma properties (see :doc:`../update_and_conv/update_and_conv`) and the
emitted spectrum (see :ref:`spectrum`).


.. _randomsampling:

Random Sampling Basics
======================

During both the initialization of these photons and their propagation through the ejecta are modeled through
probabilistic processes. This involves assigning probabilities to the occurrence of certain events or properties.
For example, during isotropic scattering, finding a photon scattering into any given direction is equally likely.
During the Monte Carlo simulation, assignments
according to these probabilities have to be frequently performed. For this purpose, so-called Random
Number Generators are available. These produce (pseudo-) random numbers
:math:`z` uniformly distributed on the interval :math:`[0,1]`. The challenge
now lies in using these numbers to sample any physical process involved in the
Radiative transfer calculation. From a probability theory point of view, this
just implies finding a mapping between the probability distribution governing the
physical process and the one underlying the Random Number Generator. This
process is typically referred to as random sampling.

Inverse transformation method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::
    This is a very superficial and sloppy description of Random sampling. More
    detailed and rigorous accounts are found in the standard literature, for
    example in :cite:`Kalos2008`.

The simplest and most-used technique in Monte Carlo radiative transfer
applications is referred to as the inverse transformation method and involves
the cumulative distribution function. In general, a random process captured by
the random variable :math:`X` is governed by a probability density
:math:`\rho_X(x)` (the continuous counterpart to discrete probabilities), with
:math:`\rho_X(x) \mathrm{d}x` describing the probability of the variable taking
values in the interval :math:`[x, x+\mathrm{d}x]`. The cumulative distribution
function in turn describes, as the name suggests, the probability of the
variable taking any value between :math:`-\infty` and :math:`x`:

.. math::

    f_X(x) = \int_{-\infty}^x \mathrm{d}x \rho_X(x)

Since the probability density is by definition always positive, the cumulative
distribution function is monotonically increasing. This constitutes the basis
for the inverse transformation function. Consider two random variables,
:math:`X` and :math:`Y`. A mapping between those may be established by equating
their cumulative distribution functions. Numbers :math:`y` distributed
according to one of the probability densities (in this case :math:`\rho_Y(y)`)
may then be used to sample the other process by

.. math::

  x = f_X^{-1}\left(f_Y(y)\right).

For the Random Number Generators described above, the cumulative distribution
function is trivial, namely :math:`f_Z(z) = z`. However, the inverse
distribution sampling method relies on finding the analytic inverse function of
the cumulative distribution function governing the physical processes to be
sampled. If this is not possible, other sampling methods, such as von-Neumann
rejection sampling techniques, have to be used.

Examples
^^^^^^^^

A few examples are provided to illustrate the random sampling process using the
inverse transformation technique.

Isotropic Scattering
--------------------

Consider the case of an isotropic scattering.
Here, all scattering angles are equally likely. Thus, the corresponding
(normalized) probability density and the cumulative distribution function are given by

.. math::

    \rho_{\mu}(\mu) &= \frac{1}{2}\\
    f_{\mu}(\mu) &= \frac{1}{2} (\mu - 1).

This leads to the sampling rule

.. math::

    \mu = 2 z - 1.

Next Interaction Location
-------------------------

The probability of a photon interacting after covering an optical depth
:math:`\tau` is given by (see :ref:`propagation`)

.. math::

    \rho_{\tau}(\tau) &= \exp(-\tau)\\
    f_{\tau}(\tau) &= 1 - \exp(-\tau).


With the inverse transformation method, the optical depth to the next interaction location may then be sampled by 

.. math::

    \tau = - \mathrm{ln}(1 - z)
  
    
which is equivalent to

.. math::

    \tau = - \mathrm{ln}z.