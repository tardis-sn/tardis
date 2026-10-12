******************
TARDIS-High Energy
******************
This code module simulates the propagation of high energy photons through the ejecta and can be used to study the transport of high energy photons within the ejecta as well as the energy deposition of high energy photons and positrons.

Thermonuclear supernovae are powered by the radioactive decay of :sup:`56` Ni and :sup:`56` Co. The photons produced by the radioactive decay often are not able to escape the ejecta until later times in the supernova when the density decreases. Once the high energy photons are able to escape they can provide information about the mass and distribution of radioactive isotopes created in the supernova.

The code can also be used to study the :math:`\gamma`-ray spectra and energy deposited by radioactive isotopes in radioactively powered transients.

For high-energy transport, the innermost radial shell is treated as a central
cell extending to the origin. An inward packet remains in that cell while it
passes through the origin in Cartesian coordinates, then crosses the shell's
outer surface on the opposite side. The origin is not treated as a boundary,
and packets are neither absorbed nor assigned a new shell there.

Interaction opacities
=====================

Detailed transport supports the ``"tardis"`` and ``"kasen"``
photoabsorption prescriptions. The Kasen prescription implements Appendix
equation A3 of :cite:`Kasen2006` using the composition moment
:math:`\sum_i n_i Z_i^5`; the default remains ``"tardis"``. Interaction types
are sampled from the comoving-frame component opacities. The Doppler factor is
applied only when converting their sum to the rest-frame path opacity.

Gamma-ray deposition outputs
============================

High-energy transport reports two complementary gamma-ray deposition
quantities. ``gamma_ray_deposited_energy`` is the event-based, cell-integrated
energy deposited in each shell and time bin, in ergs.

``gamma_ray_deposition_estimator`` is the lower-variance path estimator from
Appendix equation A10 of :cite:`Kasen2006`. For every packet segment it uses

.. math::

   \chi_{\mathrm{dep}} = \chi_{\mathrm{Compton}} F(E)
       + \chi_{\mathrm{photo}},

where :math:`F(E)` is the normalized Klein--Nishina-weighted mean fraction of
photon energy transferred to electrons, evaluated by direct angular
Gauss--Legendre quadrature of the transport's Klein--Nishina differential
cross section. In particular, :math:`F(1\,\mathrm{MeV}) \simeq 0.44`.
Segment tallies use
comoving packet energy and comoving path length. TARDIS divides the
cell-integrated tally by the time-bin width and the shell volume at the
representative time, so the returned estimator is a volumetric deposition
rate in
:math:`\mathrm{erg\,s^{-1}\,cm^{-3}}`.

The estimator includes Compton energy loss and photoabsorption as
defined in the appendix. Pair creation still contributes to the event-based
deposition output but is intentionally excluded from the estimator to match the
published equation.




.. toctree::
    packetinitialization
    opacities
    decayenergy
