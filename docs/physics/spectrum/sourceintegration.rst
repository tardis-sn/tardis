.. _formal_integral:

*****************************************
Direct integration of the radiation field
*****************************************

:cite:`Lucy1999a` describes an alternative method for the generation of
synthetic supernova spectra. Instead of using the frequency and energy of
virtual Monte Carlo (MC) packets to create a spectrum through binning, one can
formally integrate the line source functions which can be calculated from
volume-based estimators collected during the MC simulation. Spectra generated
using this method are virtually noise-free. However, statistical fluctuations
inherent to Monte Carlo calculations still affect the determination of the
source functions and thus indirectly introduce an uncertainty in the spectral
synthesis process.

.. warning::

  The current implementation of the formal integral has several limitations.
  Please consult the corresponding section below to ensure that these
  limitations do not apply to your application.


Estimators
==========

The procedure relies on determining line absorption rates, which are calculated
during the Monte Carlo simulation by incrementing the volume-based estimator

.. math::

   \dot E_{lu} = \frac{1}{\Delta t V} \left( 1- e^{-\tau_{lu}}\right) \sum
   \epsilon

Here, the sum involves all packages in a given shell that come into resonance
with the transition :math:`u \rightarrow l`, :math:`\epsilon` denotes the
energy of one such packet, and :math:`\tau_{lu}` the Sobolev optical depth of
the line transition.

After the Monte Carlo radiative transfer step, a level absorption estimator is
calculated by summing up all absorption rates for transitions which lead to the
target level

.. math::

   \dot E_u = \sum_{i < u} \dot E_{iu}

Finally, the line source function for each transition can be derived with 

.. math::

   \left( 1- e^{-\tau_{lu}}\right) S_{ul} = \frac{\lambda_{ul} t}{4 \pi} q_{ul}
   \dot E_u

Here, :math:`\lambda_{ul}` is the wavelength of the line  :math:`u \rightarrow
l`, and :math:`q_{ul}` is the corresponding branching ratio, i.e. the fraction
of de-excitations of level :math:`u` proceeding via the transition
:math:`u\rightarrow l`. For convenience, the attenuating factor is kept on the
left-hand side because it is the product of the two that will appear in later
formulae.

Finally, if the contribution by electron-scattering has to be taken into
account, estimators for the diffuse radiation field in the blue and red wings
of the different transitions are needed. These can again be determined with
volume-based estimators according to

.. math::

    J_{lu}^b = \frac{\lambda_{lu}}{4 \pi \Delta t V}  \sum \varepsilon


and

.. math::

    J_{lu}^r = J_{lu}^b e^{-\tau_{lu}} + S_{ul} (1 + e^{-\tau_{lu}})


Integration
===========

Calculating the emergent spectrum proceeds by casting rays parallel to the line
of sight to the observer through the ejecta. The distance along these rays will
be measured by :math:`z` and the offset to ejecta centre by the impact
parameter :math:`p`. The following figure illustrates this "impact geometry":

.. image:: https://i.imgur.com/WwVHp5c.png

The emergent monochromatic luminosity can then be obtained by integrating over
the limiting specific intensities of these rays

.. math::

   L_\nu  = 8 \pi^2 \int_0^\infty I_\nu (p) p dp

To obtain the limiting intensity, we have to consider the different line
resonances along the ray and calculate how much radiation is added and removed.
At the resonance point with the :math:`k`-th line transition, the intensity
increment is

.. math::

   I_k^r = I_k^b e^{-\tau_k} + \left( 1- e^{-\tau_k}\right) S_{k}

In the absence of continuum interactions, the relation

.. math::

   I_{k+1}^b = I_k^r

establishes the connection to the next resonance point. If electron-scattering
is taken into account its contribution between successive resonances has to be
considered

.. math::

   I_{k+1}^b = I_k^r + \Delta \tau_k \left[ \frac 1 2(J_k^r + J_{k+1}^b) -
   I_k^r  \right]


Thus, by recursively applying the above relations for all resonances occurring
on the ray, the limiting specific intensity for the final integration can be
calculated. The boundary conditions for this process are either :math:`I_0^r =
B_\nu(T)` if the ray intersects the photosphere or :math:`I_0^r = 0` otherwise.

Implementation Details
======================

We seek to integrate all emissions at a certain wavelength :math:`\nu` along a
ray with impact parameter :math:`p`. Because the supernova ejecta is expanding
homologously, the co-moving frame frequency is continuously shifted to longer
wavelengths due to the relativistic Doppler effect as the packet/photon
propagates.

To find out which lines can shift into the target frequency, we need to
calculate the maximum Doppler shift along a given ray. At any point, the
Doppler effect in our coordinates is

.. math::

   \nu = \nu_0 \left( 1 + \beta \mu \right)

where :math:`\beta = \frac v c`, and :math:`\mu = \cos \theta`. Here
:math:`\theta` is the angle between the radial direction and the ray to the
observer, as shown in the figure below. Because we are in the homologous
expansion regime :math:`v = \frac r t`. Solving for :math:`\nu_0` in the above
gives the relation we seek, but we require an expression for :math:`\mu`.
Examining the figure, we see that for positive :math:`z` the angle
:math:`\theta_2` can be related to the :math:`z` coordinate of the point C by

.. math::

   \cos \theta_2 = \frac{z_c}{r} = \mu


and in turn :math:`z_c` can be given as

.. math::

   z_c = \sqrt{r_c^2 + p_c^2}

where the subscripts indicate the value at point C. By symmetry, the
intersection point for negative :math:`z` is simply :math:`-z_c`.

Using the expression for :math:`\mu`, :math:`\beta` above leads to the
dependence on :math:`r` cancelling, and solving for :math:`\nu_0` gives

.. math::

   \nu_0 = \frac{\nu}{1 + \frac{z}{ct}}

For any given shell and impact parameter, we can thus find the maximum and
minimum co-moving frequency that will give the specified lab frame frequency.
This allows us to find the section of the line list with the transitions whose
resonances have to be considered in the calculation of the limiting specific
intensity.

Current Limitations
===================

The current implementation of the formal integral has some limitations:

* once electron scattering is included, the scheme only produces accurate
  results when multiple resonances occur on the rays. This is simply because
  otherwise the :math:`J^b` and :math:`J^r` do not provide an accurate
  representation of the diffuse radiation field at the current location on the
  ray. Also, :math:`d\tau` can become large which can create unphysical,
  negative intensities.

It is always advised to check the results of the formal integration against the
spectrum constructed from the emerging Monte Carlo packets.
