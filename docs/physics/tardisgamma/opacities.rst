******************
Sources of Opacity
******************
Within the ejecta there are several sources of opacity for the :math:`\gamma`-rays that we need to take into account. 

Pair Creation
=============
This form of opacity is dominant when the energy of the :math:`\gamma`-ray is greater than or equal to :math:`2m_e c^2`.
It occurs when a :math:`\gamma`-ray passes a nucleus and creates an electron-positron pair. The positron forms a bound-state with an electron called positronium which decay into two photons at 0.511 MeV or three photons. 
The photon packet keeps its energy but is sent into a new direction.
The pair-production coefficients are:

.. math::

   \alpha_{pp}(1.5 < h\nu < 1.022) = \rho [ \frac{Z_{\text{Si}}^2}{m_{\text{Si}}} (1 - X_{\text{IGE}}) + \frac{Z_{\text{Fe}}^2}{m_{\text{Fe}}} X_{\text{IGE}}]1.0063(h\nu - 1.022) \times 10^{-27}

    \alpha_{pp}(h\nu \geq 1.5) = \rho [ \frac{Z_{\text{Si}}^2}{m_{\text{Si}}} (1 - X_{\text{IGE}}) + \frac{Z_{\text{Fe}}^2}{m_{\text{Fe}}} X_{\text{IGE}}] [0.0481 + 0.301(h\nu - 1.5) \times 10^{-27}


Compton Scattering
==================
Compton scattering is inelastic scattering that occurs when a high frequency photon collides with an electron. Unlike Thomson scattering, which is elastic, during the collision some of the energy from the photon is transferred to the electron and the scattered photon is lower energy than the incident photon. 
This means that the scattered photon also has a lower frequency and a longer wavelength. In the simulation we assume that the electrons are moving significantly slower than the photons.
To find the energy of the scattered photon we use the equation we multiply the initial energy with the compton fraction :math:`f_C`

.. math::

   f_C = \frac{1}{1 + x(1 - \cos{\theta})}

Where :math:`\theta` is the angle the photon is scattered at and :math:`x = \frac{h \nu}{m_e c^2}`

The partial cross section from 0 to :math:`\theta` in terms of f is:

.. math::

   \sigma(f) = \frac{3 \sigma_T}{8x} \frac{\left( x^2 - 2x - 2 \right) \ln(f)}{x^2} + \frac{f^2 - 1}{2f^2} +  \frac{f - 1}{x} [\frac{1}{x} + 2f + \frac{1}{xf}]
   
Where f can range from 1 to 1+2x and :math:`\sigma_T` is the Thomson scattering cross section.

The scattering angle is given by:

.. math::

   \theta_C = \cos^{-1}({1 - \frac{f-1}{x}})


The total integrated Compton scattering coeffcient is:

.. math::

   \alpha_C = n_e \frac{3}{4} \sigma_T [\frac{1+x}{x^3} \frac{2x(1 + x)}{1 + 2x} - \ln(1 + 2x) + \frac{1}{2x} \ln(1 + 2x) - \frac{1 + 3x}{(1 + 2x)^2}]

The direction vector is then rotated by :math:`\theta` to get the new direction and the frequency is updated.

Photoabsorption
===============
This occurs when the photon is completely absorbed by a material.
The coefficient is:

.. math::

   \alpha_{pa}(\nu) = 1.16 \times 10^{-24} (h\nu)^{-3.13}  \frac{\rho}{m_{\text{Si}}} (1 - X_{\text{IGE}}) + 25.7 \times 10^{-24} (h\nu)^{-3}  \frac{\rho}{m_{\text{Fe}}} X_{\text{IGE}}
