LTE Plasma
----------

The `LTEPlasma` plasma class is the child of `BasePlasma` but is the first class that actually calculates plasma conditions.
After running exactly through the same steps as `BasePlasma`, `LTEPlasma` will start calculating the `partition functions <http://en.wikipedia.org/wiki/Partition_function_(statistical_mechanics)>`_.

.. math::
    Z_{i, j} = \sum_{k=0}^{max (k)} g_k \times e^{-E_k / (k_\textrm{b} T)}

where Z is the partition function, g is the degeneracy factor, E the energy of the level and T the temperature of the radiation field.

The next step is to calculate the ionization balance using the `Saha ionization equation <http://en.wikipedia.org/wiki/Saha_ionization_equation>`_.
and then calculating the Number density of the ions (and an electron number density) in a second step.
First :math:`g_e=\left(\frac{2 \pi m_e k_\textrm{B}T_\textrm{rad}}{h^2}\right)^{3/2}` is calculated (in `LTEPlasma.update_t_rad`),
followed by calculating the ion fractions (`LTEPlasma.calculate_saha`).

.. math::

    \frac{N_{i, j+1}\times N_e}{N_{i, j}} &= \Phi_{i, (j+1)/j} \\
    \Phi_{i, (j+1)/j} &= g_e \times \frac{Z_{i, j+1}}{Z_{i, j}} e^{-\chi_{j\rightarrow j+1}/k_\textrm{B}T}\\

In the second step (`LTEPlasma.calculate_ionization_balance`), we calculate in an iterative process the electron density and the number density for each ion species.

.. math::
    N(X) &= N_1 + N_2 + N_3 + \dots\\
    N(X) &= N_1 + \frac{N_2}{N_1} N_1 + \frac{N_3}{N_2}\frac{N_2}{N_1} N_1 + \frac{N_4}{N_3}\frac{N_3}{N_2}\frac{N_2}{N_1} N_1 + \dots\\
    N(X) &= N_1 (1 + \frac{N_2}{N_1} + \frac{N_3}{N_2}\frac{N_2}{N_1} + \frac{N_4}{N_3}\frac{N_3}{N_2}\frac{N_2}{N_1} + \dots)\\
    N(X) &= N_1 \underbrace{(1 + \frac{\Phi_{i, 2/1}}{N_e} + \frac{\Phi_{i, 2/2}}{N_e}\frac{\Phi_{i, 2/1}}{N_e} +
            \frac{\Phi_{i, 4/3}}{N_e}\frac{\Phi_{i, 3/2}}{N_e}\frac{\Phi_{i, 2/1}}{N_e} + \dots)}_{\alpha}\\
    N_1 &= \frac{N(X)}{\alpha}

Initially, we set the electron density (:math:`N_e`) to the sum of all atom number densities. After having calculated the
ion species number densities, we recalculate the electron density by weighting the ion species number densities with their
ion number (e.g. neutral ion number densities don't contribute at all to the electron number density, once ionized contribute with a
factor of 1, twice ionized contribute with a factor of two, ....).

Finally, we calculate the level populations (`LTEPlasma.calculate_level_populations`) by using the calculated ion species number densities:

.. math::
    N_{i, j, k} = \frac{g_k}{Z_{i, j}}\times N_{i, j} \times e^{-\beta_\textrm{rad} E_k}

This concludes the calculation of the plasma. In the code, the next step is calculating the :math:`\tau_\textrm{Sobolev}` using
the quantities calculated here.

Example Calculations
^^^^^^^^^^^^^^^^^^^^

.. .. plot:: physics/pyplot/lte_ionization_balance.py
    :include-source:

