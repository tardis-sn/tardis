***********************************
Volume-based Monte Carlo Estimators
***********************************

An integral part of Monte Carlo radiative transfer calculations consists of
reconstructing macroscopic radiation field properties from the ensemble of
packet interaction histories. TARDIS uses volume-based Monte Carlo estimators
for this purpose. This concept was originally developed by :cite:`Lucy1999` and
successively refined, for example, by :cite:`Lucy1999a`, :cite:`Lucy2002` and
:cite:`Lucy2003`.

Basic Principle
===============

The basic principle underlying volume-based estimators is best illustrated
by the example of reconstructing the mean radiation field energy density within
a certain control volume (in most cases, this will be a grid cell). A simple
approach would involve explicitly counting the number of Monte Carlo packets
which reside in the control volume at a certain time. Although being intuitive,
this approach suffers heavily from Monte Carlo shot noise, in particular in
situations in which only few (or even zero) packets reside in the control
volume when performing the counting.

The volume-based estimator approach addresses this weakness. It considers a
certain time interval :math:`\Delta t` and accounts for all packets which at
any point during the control interval resided in the control volume. Each
packet then contributes to the total time and volume-averaged radiation field
energy density according to its energy weighted by the time it spent in the
control volume :math:`\delta t` relative to the control interval :math:`\Delta
t`:

.. math::

    E = \frac{1}{\Delta V} \sum_i \varepsilon_i \frac{\delta t_i}{\Delta t}

This may be turned into the following estimator involving the trajectory
segment lengths :math:`l_i`

.. math::

    E = \frac{1}{c \Delta V \Delta t} \sum_i \varepsilon_i l_i

by exploiting :math:`l_i = \delta t_i c`. Compared to the simple counting
approach, the volume-based estimator will yield non-zero results as long as at
least one packet passed through the control volume. Additionally, as one packet
may generally contribute multiple times to the estimator (e.g. if it is
deflected by a scattering), these estimators typically suffer from less Monte
Carlo noise than simple counting reconstruction schemes.

Volume-based Estimators in TARDIS
=================================

Within TARDIS, volume-based estimators are used to reconstruct the mean
intensity and the mean frequency of the radiation field (see also
:cite:`Kerzendorf2014`). An estimator for the former is easily formulated by
using the fundamental relation :math:`J = \frac{c}{4\pi} E` and the energy density
estimator derived above:

.. math::

    J = \frac{1}{4\pi \Delta V \Delta t}\sum_i \varepsilon_i l_i D_{\mu}

An intensity-weighted estimate for the mean frequency is obtained from

.. math::

    \bar \nu = \frac{1}{4\pi \Delta V \Delta t}\sum_i \varepsilon_i \nu_i l_i D_{\mu}.

.. note::

    Compared to the estimators derived in the previous section, the ones
    presented here involve a relativistic factor, :math:`D_{\mu} = (1 - \beta
    \mu)`, which ensures the correct frame transformation behaviour of the
    estimators (to first order in :math:`v/c`).


Using the estimators just derived the radiation temperature (which should be
interpreted as a colour temperature) of the radiation field,

.. math::
    
    T_{\mathrm{R}} = \frac{h}{k_{\mathrm{B}}} \frac{\pi^4}{360 \zeta(5)} \frac{\bar \nu}{J}

may be derived. The normalization constants, involving Riemann's zeta function,
are a consequence from the definition of this colour temperature: This should
be the temperature of a black-body radiation field whose mean frequency is
equal to the one reconstructed from the Monte Carlo simulation. With the
temperature determined, the dilution factor, describing the deviation of the
radiation field from a thermal field with the same colour temperature may be calculated

.. math::

    W = \frac{\pi J}{\sigma_{\mathrm{R}} T_{\mathrm{R}}^4}

    
These two quantities, :math:`T_{\mathrm{R}}` and :math:`W` are vital for the
calculation of the plasma state of the ejecta (see :ref:`plasma_calculations`).

Finally, TARDIS also reconstructs the mean intensity of the radiation field in
the blue wing of each line transition :math:`l \rightarrow u`, which is used in
the Macro Atom treatment and in the ionisation calculation.

.. math::

    J_{lu}^b = \frac{1}{4\pi \Delta V \Delta t} \frac{t_{\mathrm{exp}}}{c} \sum_i \frac{\varepsilon_i}{\nu_{lu}} D_{\mu}.
    
The summation here only involves packets which passed through the Sobolev point
(see :ref:`Propagation <propagation>`) of the transition. For a derivation of this
estimator, in particular, for a motivation of the expansion factor involving
the time since explosion :math:`t_{\mathrm{exp}}`, we refer to
:cite:`Lucy2003`, section 6.2.
