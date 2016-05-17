******************
Packet Propagation
******************

The bulk of a Monte Carlo Radiative Transfer calculation is spent on
determining the propagation history of the different packets. After a packet is
initialised, it is launched and may then perform interactions with the
surrounding material. This occurs again in a probabilistic manner. The packet
propagation is followed until it escapes through the outer boundary of the
computational domain at which point the packet contributes to the synthetic
spectrum, the main product of a TARDIS calculation. The different spectral
features are simply a combined product of the changes in the packet properties
induced in the radiation-matter interactions.

Initialization
==============

During each TARDIS Monte Carlo simulation cycle, a large number :math:`N` of
Monte Carlo packets is initialised at the lower boundary of the computational
(i.e.  the photosphere). Since the inner boundary is currently treated as a
black-body in TARDIS, :math:`N` packets with energies

.. math::
    \varepsilon = \frac{4 R_{\mathrm{phot}}^2 \sigma_{\mathrm{R}} T_{\mathrm{phot}}^4}{N}

are initialised (the black body temperature :math:`T_{\mathrm{phot}}`, the
photospheric radius :math:`R_{\mathrm{phot}}` and the Stefan-Boltzmann constant
:math:`\sigma_{\mathrm{R}}` appear here). To commence the packet propagation,
each packet is assigned a initial propagation direction (recall that :math:`\mu
= \cos \theta` with `\theta` being the angle enclosed by the photon path and
the radial direction)

.. math::
    \mu = \sqrt{z}

and an initial frequency :math:`\nu` in random number experiments, using a
Random number generator which provides uniformly distributed numbers :math:`z`
on the interval :math:`[0,1]`. The frequency assignment is more involved than
selecting an initial propagation direction, since the Planck function has to be
sampled. TARDIS used the technique described in :cite:`Carter1975` and
summarized in :cite:`Bjorkmann1999` for this purpose.

Propagation in a spherical domain
=================================

Once the initial packet properties are assigned, the propagation process
commences. Without interacting, a packet, like a photon, will propagate on a
straight trajectory.

.. note::
    Since TARDIS is designed for systems for which the Newtonian limit may be
    safely applied, no general relatistic effects which would force photons on
    curved trajectories are included.

In spherical geometry, this propagation process may be illustrated by the
following sketch (taken from :cite:`Noebauer2014`):


.. image::
    ../graphics/spherical_symmetry.png
    :width: 400


The packet starts the propagation at :math:`r_i` along the direction
:math:`\mu_i`. After covering a distance :math:`l`, the packet is now located
at

.. math::
    r_f = \sqrt{r_i^2 + l^2 + 2 l r_i \mu_i}.

Note that the propagation direction has also changed and now takes the value 

.. math::
    \mu_f = \frac{l + r_i \mu_i}{r_f}.

Numerical and Physical Events
=============================

As a packet propagates through the computational domain, a number of events may
trigger changes in the packet properties. Naturally, physical radiation-matter
interactions are such events. These always occur after the packet has covered a
distance corresponding to the optical depth

.. math::

    \tau = -\log z,

which is again assigned probabilistically, in accordance with the stochastic
nature of the Monte Carlo approach. Translating the optical depth to the next
interaction into a physical distance is not straight-forward in the presence of
frequency-dependent interaction process such as atomic line interactions. The
detailed procedure is outlined in the following section.

In addition to the physical processes, numerical events which are a consequence
of the spatial discretization of the computational domain require interrupting
the propagating process. In TARDIS, as in many other numerical codes, physical
quantities are stored on a discrete mesh. Within the different cells, which in
our case are radial shells, these discrete values determine a (spatially)
constant plasma state. As a consequence, whenever a packet propagates into a
new cell, important quantities which are relevant for performing
radiation-matter interactions have to be re-evaluated in accordance with the
new state of the ambient material. Thus, during the packet propagation, the
distance to the next radial shell is tracked to predict when the packet crosses
into a new shell. Special care is taken at the edges of the computational
domain. If a packet crosses back into the photosphere, it is discarded. Its
propagation is stopped and it is no longer considered. Instead processing the
next packet of the population is started. Similarly, the propagation is stopped
if the packet escapes through the outer surface of the domain. However, in this
case the packet contributes to the final emergent spectrum (see :doc:`Spectrum
Formation <virtualpackets>`).

Physical Events
===============

As noted above, translating the optical depth, which determines when the next
physical interaction occurs, is non-trivial as soon as frequency-dependent
processes are considered. Currently, TARDIS incorporates the electron
scatterings and interactions with atomic line transitions. These two
interactions mechanisms constitute the main sources of opacity in Type Ia
supernovae. Since the main focus of TARDIS is to calculate optical spectra,
electron-scatterings are treated in the elastic low-energy limit as classical
Thomson scatterings.


- optical depth summation

- Thomson scattering

- Resonant Line Interaction


Implementation: Main Propagation loop
=====================================

In summary of the concepts outlined above, the main Monte Carlo process within
TARDIS consists of successively processing all packets with represent the
radiation field emitted by the photosphere in the following way:

* initialize the packet: assign initial energy, direction and frequency
* launch the packet: now the propagation of this packet is followed until one of the termination events is triggered
* follow the propagation:
    * calculate the distance to the next shell and determine the distance to the next physical interaction
    * the packet covers the shorter of these two distances:
       * if the new shell is reached first, propagate into the shell and recalculate both distances
       * if the packet has crossed through the inner domain boundary (photosphere), terminate propagation
       * likewise in case the packet escapes through outer boundary (ejecta surface): account for contribution to emergent spectrum and terminate propagation
       * if the interaction location is reached first, propagated to this location, perform interaction and recalculate both distances
    * repeat this procedure until one of the two termination events occurs

The following flow chart summarizes this process again:

*Coming Soon*
