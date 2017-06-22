******************************** Direct integration of the radiation field
********************************

.. note::

  The current implementation only works with the downbranch line interaction
  scheme.


:cite:`Lucy1999a` describes an alternative method for the generation of
synthetic supernova spectra. Instead of using the frequency and energy of
virtual Monte Carlo (MC) packets to create a spectrum through binning, one can
formally integrate the line source functions which can be calculated from
volume based estimators collected during the MC simulation. Spectra
generated using this method do not contain MC noise directly. Here the
MC nature of the simulation only affects the strengths of lines and
thus the spectra appear to be of better quality.

The procedure uses a line absorption rate estimator that is collected during
the MC simulation:

.. math::

   \dot E_{lu} = \frac{1}{\Delta t V} \left( 1- e^{-\tau_{lu}}\right) \sum
   \epsilon

where the sum is over all the packages in a given shell that come into
resonance with the transition :math:`u \rightarrow l` during the MC
run, :math:`\epsilon` is the energy of one such packet, and :math:`\tau_{lu}`
the Sobolev optical depth.

After the final MC step, a level absorption estimator is calculated,
which includes all levels which lie below the target level:

.. math::

   \dot E_u = \sum_{i < u} \dot E_{iu}

The source function for each line can then be derived from the relation

.. math::

   \left( 1- e^{-\tau_{lu}}\right) S_{ul} = \frac{\lambda_{ul} t}{4 \pi} q_{ul}
   \dot E_u

where :math:`\lambda_{ul}` is the wavelength of each line  :math:`u \rightarrow
l`, and :math:`q_{ul}` is the corresponding branching ratio. The attenuating
factor is kept on the left hand side because it is the product of the two that
will appear in later formulae.

The formal integration is based on the so-called
"elementary supernova" model, which is described in detail in Jeffery & Branch
1990. The final integral is given as

.. math::

   L_\nu  = 8 \pi^2 \int_0^\infty I_\nu (p) p dp

where :math:`p` is the impact parameter of a ray trough the supernova envelope
that reaches the distant observer, and :math:`I_\nu (p)` is the intensity along
one such ray, given by recursing through the list of attenuated source functions
from the blue to the red and adding up contributions. The relation linking the
intensity before the k:th transition :math:`u \rightarrow l` to the intensity
after is

.. math::

   I_k^r = I_k^b e^{-\tau_k} + \left( 1- e^{-\tau_k}\right) S_{k}

where the superscripts are crucial, with :math:`r` and :math:`b` referencing
the red and blue sides of the k:th transition respectively. To go from the red
side of a line to the blue side of the next we can either ignore continuum
sources of opacity, in which case

.. math::

   I_{k+1}^b = I_k^r

.. note::

   Currently the code does not perform the steps necessary to include continuum
   sources of opacity.

or include them, then requiring we perform

.. math::

   I_{k+1}^b = I_k^r + \Delta \tau_k \left[ \frac 1 2(J_k^r + J_{k+1}^b) -
   I_k^r  \right]

The starting condition for the blue to red side transition is either
:math:`I_0^r = B_\nu(T)` if the ray intersects the photosphere and :math:`I_0^r
= 0` otherwise.

We seek to integrate all emissions at a certain wavelength :math:`\nu` along a
ray with impact parameter :math:`p`. Because the supernova ejecta is expanding
homologously, the co-moving frame frequency is continuously shifted to longer
wavelength due to the relativistic Doppler effect as the packet/photon
propagates.


To find out, which lines can shift into the target frequency, we need to calculate
the maximum Doppler shift along a given ray. At any point, the Doppler effect
in our coordinates is

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

.. image:: https://i.imgur.com/WwVHp5c.png

and in turn :math:`z_c` can be given as

.. math::

   z_c = \sqrt{r_c^2 + p_c^2}

where the subscripts indicate the value at point C. By symmetry the
intersection point for negative :math:`z` is simply :math:`-z_c`.

Using the expression for :math:`\mu`, :math:`\beta` above leads to the
dependence on :math:`r` cancelling, and solving for :math:`\nu_0` gives

.. math::

   \nu_0 = \frac{\nu}{1 + \frac{z}{ct}}

For any given shell and impact parameter we can thus find the maximum and
minimum co-moving frequency that will give the specified lab frame frequency.
