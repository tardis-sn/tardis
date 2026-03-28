.. _nebular_plasma:

Nebular Plasma
--------------

The `NebularPlasma` class is a more complex description of the Plasma state than the `LTEPlasma`. It takes a dilution factor
(W) into account, which deals with the dilution of the radiation field due to geometric, line-blocking and other effects.


The calculations follow the same steps as `LTEPlasma`; however, the calculations are different and often take into account
if a particular level is :term:`meta-stable` or not.
`NebularPlasma` will start calculating the `partition functions <http://en.wikipedia.org/wiki/Partition_function_(statistical_mechanics)>`_.

.. math::

    Z_{i,j} = \underbrace{\sum_{k=0}^{max(k)_{i,j}} g_k \times e^{-E_k / (k_\textrm{b} T)}}_\textrm{metastable levels} +
            \underbrace{W\times\sum_{k=0}^{max(k)_{i,j}} g_k \times e^{-E_k / (k_\textrm{b} T)}}_\textrm{non-metastable levels}

where Z is the partition function, g is the degeneracy factor, E the energy of the level, T the temperature of the radiation field
and W the dilution factor.

The next step is to calculate the ionization balance using the `Saha ionization equation <http://en.wikipedia.org/wiki/Saha_ionization_equation>`_
and then calculate the number density of the ions (and an electron number density) in a second step.
In the first step, we calculate the ionization balance using the LTE approximation (:math:`\Phi_{i, j}(\textrm{LTE})`). Then we adjust the ionization balance using
two factors: :math:`\zeta` and :math:`\delta`.


.. _calc_zeta_label:

Calculating Zeta
^^^^^^^^^^^^^^^^

:math:`\zeta` is read in for specific temperatures and then interpolated for the target temperature.


Calculating Delta
^^^^^^^^^^^^^^^^^

:math:`\delta` is a radiation field correction factors which is calculated according to Mazzali & Lucy 1993 (:cite:`Mazzali1993`; henceforth ML93)

In ML93, the radiation field correction factor is denoted as :math:`\delta` and is calculated in Formula 15 & 20.

The radiation correction factor changes according to a ionization energy threshold :math:`\chi_\textrm{T}`
and the species ionization threshold (from the ground state) :math:`\chi_0`.

**For** :math:`\chi_\textrm{T} \ge \chi_0`

.. math::
    \delta = \frac{T_\textrm{e}}{b_1 W T_\textrm{R}} \exp(\frac{\chi_\textrm{T}}{k T_\textrm{R}} -
    \frac{\chi_0}{k T_\textrm{e}})

**For** :math:`\chi_\textrm{T} < \chi_0`

.. math::
    \delta = 1 - \exp(\frac{\chi_\textrm{T}}{k T_\textrm{R}} - \frac{\chi_0}{k T_\textrm{R}}) +
    \frac{T_\textrm{e}}{b_1 W T_\textrm{R}} \exp(\frac{\chi_\textrm{T}}{k T_\textrm{R}} -
    \frac{\chi_0}{k T_\textrm{e}}),

where :math:`T_\textrm{R}` is the radiation field temperature, :math:`T_\textrm{e}` is the electron temperature and W is the
dilution factor.


Now, we can calculate the ionization balance using equation 14 in :cite:`Mazzali1993`:

.. math::
    \Phi_{i,j} &= \frac{N_{i, j+1} n_e}{N_{i, j}} \\

    \Phi_{i, j} &= W \times[\delta \zeta + W ( 1 - \zeta)] \left(\frac{T_\textrm{e}}{T_\textrm{R}}\right)^{1/2}
    \Phi_{i, j}(\textrm{LTE}) \\


In the last step, we calculate the ion number densities according using the methods in :class:`LTEPlasma`

Finally, we calculate the level populations (:func:`NebularPlasma.calculate_level_populations`),
by using the calculated ion species number densities:

.. math::

    N_{i, j, k}(\textrm{not metastable}) &= W\frac{g_k}{Z_{i, j}}\times N_{i, j} \times e^{-\beta_\textrm{rad} E_k} \\
    N_{i, j, k}(\textrm{metastable}) &= \frac{g_k}{Z_{i, j}}\times N_{i, j} \times e^{-\beta_\textrm{rad} E_k} \\


This concludes the calculation of the nebular plasma. In the code, the next step is calculating the :math:`\tau_\textrm{Sobolev}` using
the quantities calculated here.

Example Calculations
^^^^^^^^^^^^^^^^^^^^

.. .. plot:: physics/pyplot/nebular_ionization_balance.py
    :include-source:

