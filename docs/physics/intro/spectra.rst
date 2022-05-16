.. _spectra:

*********************
Radiation and Spectra
*********************

As mentioned, one of the main tasks of TARDIS is determining the spectrum of light emitted by a supernova. In this section, we explore important quantities related to spectra.

.. note::
    
    This page is intended as a basic overview. For more detail see ???.


Important Quantities
====================

Luminosity
----------

The first quantity is related to how much light *in total* is emitted by some object -- i.e. across the entire surface over every frequency. This is called the **luminosity** :math:`L`, and is defined as the rate at which light energy is emitted by an object. So, if the total energy of light emitted in a time :math:`\Delta t` is :math:`E`, we have

.. math:: L = \frac{E}{\Delta t}.

Luminosity Density
------------------

When we are looking at a spectrum, we are not typically as interested in the total luminosity as we are interested in the luminosity at particular wavelengths/frequencies of light. However, this presents a problem: the possible frequencies of light is a continuum, while there is a finite amount of light that comes out of a source. So most frequencies won't have any photons with the *exact* same frequency, and there is a nearly zero probability that any two photons will have the *exact* same frequency (i.e. one could have a frequency of :math:`10^{14}` Hz and the other :math:`10^14+0.001` Hz, for example). This problem is demonstrated here??? using an example from TARDIS (note that in TARDIS we use finite packets of light as opposed to photons, but the concept is the same).

The best we can do is consider the luminosity within a certain range of frequencies. This inspires the concept of luminosity density, also called specific luminosity, denoted :math:`L_\nu`. This is the luminosity **per unit frequency**. So, for example, if there is a luminosity of :math:`4\times 10^{40}` erg/s between the frequencies of :math:`10^{14}` Hz and the other :math:`10^14+2` Hz, the luminosity density will be :math:`\frac{4\times 10^{40} \mathrm{erg/s}}{2\mathrm{Hz}}=2\times 10^{40}\mathrm{erg/s/Hz}` since the luminosity is within a frequency range of 2 Hz. This definition makes it so, for example, a luminosity of :math:`4\times 10^{40}` erg/s within a range of 2 Hz is the same as a luminosity of :math:`2\times 10^{40}` erg/s within a range of 1 Hz, as we might expect. That is, the luminosity density does not depend on how large or small of frequency intervals we consider, though we may lose precision if we make our intervals too small.

The luminosity density :math:`L_\nu` at a frequency :math:`\nu` is the luminosity density in a small interval containing the :math:`\nu`. We then write that the luminosity :math:`dL` between frequencies :math:`\nu` and :math:`\nu+d\nu` is

.. math:: dL=L_\nu d_\nu.

To get the total luminosity, we add up (integrate) the luminosity overr every interval, so

.. math:: L = \int L_\nu d\nu.

We have just defined luminosity density with respect to frequency. However, we can also define luminosity density with respect to wavelength :math:`L_\lambda`, which is luminosity per unit wavelength. This is done the exact same way, but by looking at the luminosity across a certain range of wavelengths of light and dividing by the width of that range. The plots of these two luminosity densities will look different due to the inverse relationship between wavelength and frequency (see ???), and both are often used to plot spectra. Note that just as before we have the total luminosity

.. math:: L = \int L_\lambda d\lambda.


.. _intensity:

Intensity
---------

Specific intensity :math:`I_\nu`, often simply refered to as intensity, is the luminosity density carried along a specific ray of light at a specific location -- technically, it is the luminosity carried along a ray per unit area per unit solid angle ???hyperlink???. This is similar to why luminosity density is luminosity per unit wavelength or frequency -- there are a continuum of possible directions and places on the surface of the object for light to be emitted.

The luminosity density :math:`dL_\nu` of a beam of light with intensity :math:`I_\nu`, with light emitted from a small area :math:`dA` over a small range of solid angles :math:`d\Omega`, will be

.. math:: dL_\nu = I_\nu d\Omega dA.

Thus if we add up the intensities of all rays of light coming from all parts of the surface going in all directions, we get the total luminosity density:

.. math:: L_\nu = \int I_\nu d\Omega dA.

.. note::

    We only show intensity and the rest of the quantities on this page using luminosity density with respect to frequency, as this is typically how they are used in TARDIS and astronomy in general???i think.

???????????? make a note about how this relates to number of photons

Specific Mean Intensity
-----------------------

Specific mean intensity (denoted by :math:`J_\nu`), is another useful quantity that tells us about the "overall intensity" of light of a frequency :math:`\nu` at a certain location, independent of the exact ray of the light. To do this, we take the average intensity at some location over all possible directions/solid angles.

Specifically, we integrate the intensity over all solid angles and divide by :math:`4\pi`, the total solid angle:

.. math:: J_\nu = \frac{\int I_\nu d\Omega}{4\pi}.

Note that now we have simply luminosity density per unit area.

Since we have :math:`\int I_\nu d\Omega= 4\pi J_\nu`, from the previous section we see the relationship between luminosity density and specific mean intensity:

.. math:: L_\nu = 4\pi\int J_\nu dA.

Integrated Mean Intensity
-------------------------

Next, we can determine the "overall intensity" of light reguardless of direction *or* frequency. This is called the integrated mean intensity, or simply the mean intensity, denoted by :math:`J`. The relationship between integrated mean intensity and specific mean intensity is the same as the relationship between the relationship between luminosity and luminosity density:

.. math:: J = \int J_\nu d\nu.


Energy Density
--------------

Our final quantity of importance is specific energy density :math:`u_\nu`, which is light energy per unit volume per unit frequency, and integrated energy density :math:`u` which is light energy per unit volume. Once again, these are related by

.. math:: u = \int u_\nu d\nu.

We interestingly have a special relationship between mean intensity and energy density: 

.. math:: u_\nu=\frac{4\pi}{c}J_\nu
    
and thus

.. math:: u=\frac{4\pi}{c}J.

This is explained, for example in ???.


Blackbody Spectra
=================

One of the major advancements in physics in the early 1900s was the discovery that all objects emit light in a way that depends on its temperature. This is called **blackbody radiation** and could be thought as the most "generic" spectrum -- it is emitted by everything and depends soley on the temperature. We will present the value of some of the above quantities for the special case of blackbody radiation from a sphere with radius :math:`R` temperature :math:`T`. For more on the derivation of these equations see ???.

To start, the luminosity density with respect to frequency :math:`\nu` is given by the Planck spectrum

.. math:: L_\nu =\frac{8\pi R^2 h\nu^3}{c^2}\frac{1}{\exp\left(\frac{h\nu}{k_BT}\right)-1}

where :math:`h` is Planck's constant, :math:`c` is the speed of light, and :math:`k_B` is Boltzmann's constant. Using the wavelength-frequency relation (stated here???) :math:`\lambda=\frac{c}{\nu}`, we can get the luminosity density with respect to wavelength:

.. math:: L_\lambda =\frac{8\pi R^2 hc}{\lambda^3}\frac{1}{\exp\left(\frac{hc}{\lambda k_BT}\right)-1}.

The Planck spectrum is plotted below for various temperatures.

.. planck spectrum image

Now, we can integrate either luminosity density to get the total luminosity of our blackbody:

.. math:: L = 4\pi R^2 \sigma_R T^4

where :math:`\sigma_R` is the Stefan-Boltzmann constant, which is obtained by doing the integral.

Another key property of blackbody radiation is that it is isotropic (the same in all directions). Every ray emitted from the blackbody thus has an equal intensity, given by the Planck function:

.. math:: I_\nu = B_\nu \equiv \frac{2 h\nu^3}{c^2}\frac{1}{\exp\left(\frac{h\nu}{k_BT}\right)-1}.


Interpreting Spectra
====================

Blackbody spectra are extremely important for understanding spectra in general. An example spectrum of ??? is shown below. It has the general shape of a blackbody -- after all, any object wth a temperature exmits a blackbody spectrum. Dips that you see at certain wavelengths??? has to do with light-matter interactions happening at that wavelength???, especially line interactions (see ???). This tells us a lot about what atoms the object emitting the spectrum is made out of, since every atom has different wavelengths??? of light that it interacts with. However, you may notice that these are not all sharp dips -- they cover a whole range of wavelengths??? even though line interactions happen at only specific wavelengths???. This is a phenominon called line broadening, specifically **doppler broadening** that causes these atomic lines to cover many frequencies. This has to do with the doppler effect caused by the atoms all moving at different velocities. Much of TARDIS deals with this doppler broadening -- the physics that is occuring is explained within the context of TARDIS here???.

.. example spectrum image