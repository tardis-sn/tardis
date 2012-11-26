Plasma
******
.. currentmodule:: tardis.plasma

This module calculates the ionization balance and level populations in the Plasma, give a abundance fraction, temperature
and density. After calculating the state of the plasma, these classes are able to calculate :math:`\tau_\textrm{sobolev}`
for the supernova radiative transfer. The simplest Plasma (`Plasma`) only calculates the atom number densities, but serves
as a base for all Plasma classes. The next more complex class is `LTEPlasma` which will calculate the aforementioned quantities in
Local Thermal Equilibrium conditions (LTE). The `NebularPlasma`-class inherits from `LTEPlasma` and uses a more complex
description of the Plasma (for details see `Nebular Plasma`_).

.. important::
    In this documentation we use the indices :math:`i, j, k` to mean atomic number, ion number and level number respectively.



Base Plasma
===========

`Plasma` serves as the base class for all plasmas and can just calculate the atom number densities for a given input of
abundance fraction.

.. math::
    N_{atom} = \rho_\textrm{total} \times \textrm{Abundance fraction} / m_\textrm{atom}

In the next step the line and level tables are purged of entries that are not represented in the
abundance fractions are saved in `Plasma.levels` and `Plasma.lines`. Finally, the function `Plasma.update_t_rad` is called
at the end of initialization to update the plasma conditions to a new :math:`T_\textrm{radiation field}` (with the give t_rad).
This function is the same in the other plasma classes and does the main part of the calculation. In the case of `Plasma` this is only
setting `Plasma.beta_rad` to :math:`\frac{1}{k_\textrm{B}T_\textrm{rad}}`.


Here's an example how to instantiate a simple base plasma::


    >>> from tardis import atomic, plasma
    >>> atom_data = atomic.AtomData.from_hdf5()
    >>> my_plasma = plasma.Plasma({'Fe':0.5, 'Ni':0.5}, 10000, 1e-13, atom_data)
    >>> print my_plasma.abundances
    atomic_number abundance_fraction number_density
    ------------- ------------------ --------------
               28                0.5  513016973.936
               26                0.5  539183641.472



LTE Plasma
==========

The `LTEPlasma` plasma class is the child of `Plasma` but is the first class that actually calculates plasma conditions.
After running exactley through the same steps as `Plasma`, `LTEPlasma` will start calculating the `partition functions <http://en.wikipedia.org/wiki/Partition_function_(statistical_mechanics)>`_.

.. math::
    Z_{i, j} = \sum_{k=0}^{max (k)} g_k \times e^{-E_k / (k_\textrm{b} T)}

, where Z is the partition function, g is the degeneracy factor, E the energy of the level and T the temperature of the radiation field.

The next step is to calculate the ionization balance using the `Saha ionization equation <http://en.wikipedia.org/wiki/Saha_ionization_equation>`_.
and then calculating the Number density of the ions (and an electron number density) in a second step.
First :math:`g_e=\left(\frac{2 \pi m_e k_\textrm{B}T_\textrm{rad}}{h^2}\right)^{3/2}` is calculated, followed by calculating
the ion fractions.

.. math::

    \frac{N_{i, j+1}\times N_e}{N_{i, j}} &= \Phi_{i, (j+1)/j} \\
    \Phi_{i, (j+1)/j} &= g_e \times \frac{Z_{i, j+1}}{Z_{i, j}} e^{-\chi_{j\rightarrow j+1}/k_\textrm{B}T}\\

In the second step, we calculate in an iterative process the electron density and the number density for each ion species.

.. math::
    N(X) &= N_1 + N_2 + N_3 + \dots\\
    N(X) &= N_1 + \frac{N_2}{N_1} N_1 + \frac{N_3}{N_2}\frac{N_2}{N_1} N_1 + \frac{N_4}{N_3}\frac{N_3}{N_2}\frac{N_2}{N_1} N_1 + \dots\\
    N(X) &= N_1 (1 + \frac{N_2}{N_1} + \frac{N_3}{N_2}\frac{N_2}{N_1} + \frac{N_4}{N_3}\frac{N_3}{N_2}\frac{N_2}{N_1} + \dots)\\
    N(X) &= N_1 \underbrace{(1 + \frac{\Phi_{i, 2/1}}{N_e} + \frac{\Phi_{i, 2/2}}{N_e}\frac{\Phi_{i, 2/1}}{N_e} +
            \frac{\Phi_{i, 4/3}}{N_e}\frac{\Phi_{i, 3/2}}{N_e}\frac{\Phi_{i, 2/1}}{N_e} + \dots)}_{\alpha}\\
    N_1 &= \frac{N(X)}{\alpha}

Initially, we set the electron density (:math:`N_e`) to the sum of all atom number densities. After having calculated the
ion species number densities we recalculated the electron density by weighting the ion species number densities with their
ion number (e.g. neutral ion number densities don't contribute at all to the electron number density, once ionized contribute with a
factor of 1, twice ionized contribute with a factor of two, ....).

Finally we calculate the level populations, by using the calculated ion species number densities:

.. math::
    N_{i, j, k} = \frac{g_k}{Z_{i, j}}\times N_{i, j} \times e^{-\beta_\textrm{rad} E_k}

This concludes the calculation of the plasma. In the code, the next step is calculating the :math:`\tau_\textrm{Sobolev}` using
the quantities calculated here.


Nebular Plasma
==============





Level Populations
-----------------

We first calculate the number densities for each atom


Next we calculate the partition function (Z) where i is the level number and j the ion number:


.. automodapi:: tardis.plasma
