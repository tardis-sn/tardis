BasePlasma
******
.. currentmodule:: tardis.plasma

This module calculates the ionization balance and level populations in the BasePlasma, give a abundance fraction, temperature
and density. After calculating the state of the plasma, these classes are able to calculate :math:`\tau_\textrm{sobolev}`
for the supernova radiative transfer. The simplest BasePlasma (`BasePlasma`) only calculates the atom number densities, but serves
as a base for all BasePlasma classes. The next more complex class is `LTEPlasma` which will calculate the aforementioned quantities in
Local Thermal Equilibrium conditions (LTE). The `NebularPlasma`-class inherits from `LTEPlasma` and uses a more complex
description of the BasePlasma (for details see `Nebular BasePlasma`_).

.. important::
    In this documentation we use the indices :math:`i, j, k` to mean atomic number, ion number and level number respectively.



Base BasePlasma
===========

`BasePlasma` serves as the base class for all plasmas and can just calculate the atom number densities for a given input of
abundance fraction.

.. math::
    N_{atom} = \rho_\textrm{total} \times \textrm{Abundance fraction} / m_\textrm{atom}

In the next step the line and level tables are purged of entries that are not represented in the
abundance fractions are saved in `BasePlasma.levels` and `BasePlasma.lines`. Finally, the function `BasePlasma.update_t_rad` is called
at the end of initialization to update the plasma conditions to a new :math:`T_\textrm{radiation field}` (with the give t_rad).
This function is the same in the other plasma classes and does the main part of the calculation. In the case of `BasePlasma` this is only
setting `BasePlasma.beta_rad` to :math:`\frac{1}{k_\textrm{B}T_\textrm{rad}}`.


Here's an example how to instantiate a simple base plasma::


    >>> from tardis import atomic, plasma
    >>> atom_data = atomic.AtomData.from_hdf5()
    >>> my_plasma = plasma.BasePlasma({'Fe':0.5, 'Ni':0.5}, 10000, 1e-13, atom_data)
    >>> print my_plasma.abundances
    atomic_number abundance_fraction number_density
    ------------- ------------------ --------------
               28                0.5  513016973.936
               26                0.5  539183641.472



LTE BasePlasma
==========

The `LTEPlasma` plasma class is the child of `BasePlasma` but is the first class that actually calculates plasma conditions.
After running exactley through the same steps as `BasePlasma`, `LTEPlasma` will start calculating the `partition functions <http://en.wikipedia.org/wiki/Partition_function_(statistical_mechanics)>`_.

.. math::
    Z_{i, j} = \sum_{k=0}^{max (k)} g_k \times e^{-E_k / (k_\textrm{b} T)}

, where Z is the partition function, g is the degeneracy factor, E the energy of the level and T the temperature of the radiation field.

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
ion species number densities we recalculated the electron density by weighting the ion species number densities with their
ion number (e.g. neutral ion number densities don't contribute at all to the electron number density, once ionized contribute with a
factor of 1, twice ionized contribute with a factor of two, ....).

Finally we calculate the level populations (`LTEPlasma.calculate_level_populations`), by using the calculated ion species number densities:

.. math::
    N_{i, j, k} = \frac{g_k}{Z_{i, j}}\times N_{i, j} \times e^{-\beta_\textrm{rad} E_k}

This concludes the calculation of the plasma. In the code, the next step is calculating the :math:`\tau_\textrm{Sobolev}` using
the quantities calculated here.

Here's an example::

    >>> from tardis import atomic, plasma
    >>> atom_data = atomic.AtomData.from_hdf5()
    >>> my_plasma = lasma.LTEPlasma({'Fe':0.5, 'Ni':0.5}, 10000, 1e-13, atom_data)
    >>> print my_plasma.partition_functions
    atomic_number ion_number partition_function
    ------------- ---------- ------------------
               26          0      59.6505567482
               26          1      66.8959778772
               26          2      28.1862532138
               26          3      6.56910482208
               26          4      24.4529907641
               26          5      26.1068211267
               26          6      18.4458984351
               26          7      8.60710657603
              ...        ...                ...
               28          1      19.3995555128
               28          2      20.1680211885
               28          3      26.6331919318
               28          4      23.6553761087
               28          6      20.9840329365
               28          7      21.7719521018
               28          8      15.7179889906
    >>> print my_plasma.levels # The last column has been calculated and changes when updating the t_rad
    atomic_number ion_number level_number ...  g  metastable   number_density
    ------------- ---------- ------------ ... --- ---------- ------------------
               26          0            0 ...   9       True  0.000169596945909
               26          0            1 ...   7       True  0.000124246414266
               26          0            2 ...   5       True  8.51442728653e-05
               26          0            3 ...   3       True  4.97509739812e-05
               26          0            4 ...   1       True  1.63704372998e-05
               26          0            5 ...  11       True  7.64985795918e-05
               26          0            6 ...   9       True  5.86784715998e-05
               26          0            7 ...   7       True  4.33893909491e-05
              ...        ...          ... ... ...        ...                ...
               28          8           19 ...   9      False 2.37668235834e-198
               28          8           20 ...   7      False 9.05523315872e-199
               28          8           21 ...   5      False 5.65796533522e-199
               28          8           22 ...   9      False 9.46379424911e-199
               28          8           23 ...  11      False 1.03092538574e-198
               28          8           24 ...   7      False 5.66496972859e-199
               28          8           25 ...   7      False 4.37244215215e-199

    >>> print my_plasma.update_t_rad(5000)
    atomic_number ion_number level_number ...  g  metastable number_density
    ------------- ---------- ------------ ... --- ---------- --------------
               26          0            0 ...   9       True  12715.9908741
               26          0            1 ...   7       True  8774.58025089
               26          0            2 ...   5       True  5768.96016405
               26          0            3 ...   3       True   3282.7558399
               26          0            4 ...   1       True  1066.29463495
               26          0            5 ...  11       True  2116.75561292
               26          0            6 ...   9       True  1522.20007689
               26          0            7 ...   7       True   1070.1033671
              ...        ...          ... ... ...        ...            ...
               28          8           19 ...   9      False            0.0
               28          8           20 ...   7      False            0.0
               28          8           21 ...   5      False            0.0
               28          8           22 ...   9      False            0.0
               28          8           23 ...  11      False            0.0
               28          8           24 ...   7      False            0.0
               28          8           25 ...   7      False            0.0




Nebular BasePlasma
==============

The `NebularPlasma` class is a more complex description of the BasePlasma state than the `LTEPlasma`. It takes a dilution factor
(W) into account, which deals with the dilution of the radiation field due to geometric, line-blocking and other effects.


The calculations follow the same steps as `LTEPlasma`, however the calculations are different and often take into account
if a particular level is :term:`meta-stable` or not.
`NebularPlasma` will start calculating the `partition functions <http://en.wikipedia.org/wiki/Partition_function_(statistical_mechanics)>`_.

.. math::

    Z_{i,j} = \underbrace{\sum_{k=0}^{max(k)_{i,j}} g_k \times e^{-E_k / (k_\textrm{b} T)}}_\textrm{metastable levels} +
            \underbrace{W\times\sum_{k=0}^{max(k)_{i,j}} g_k \times e^{-E_k / (k_\textrm{b} T)}}_\textrm{non-metastable levels}

, where Z is the partition function, g is the degeneracy factor, E the energy of the level, T the temperature of the radiation field
and W the dilution factor.

The next step is to calculate the ionization balance using the `Saha ionization equation <http://en.wikipedia.org/wiki/Saha_ionization_equation>`_.
and then calculating the Number density of the ions (and an electron number density) in a second step.
In the first step, we calculate the ionization balance using the LTE approximation (:math:`\Phi_{i, j}(\textrm{LTE})`). Then we adjust the ionization balance using
two factors :math:`\zeta` and :math:`\delta`.

Calculating Zeta
----------------

:math:`\zeta` is read in for specific temperatures and then interpolated for the target temperature.


Calculating Delta
-----------------

:math:`\delta` is a radiation field correction factors which is calculated according to Mazzali & Lucy 1993 (:cite:`1993A&A...279..447M`; henceforth ML93)

In ML93 the radiation field correction factor is denoted as :math:`\delta` and is calculated in Formula 15 & 20

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

where :math:`T_\textrm{R}` is the radiation field Temperature, :math:`T_\textrm{e}` is the electron temperature and W is the
dilution factor.


Now we can calculate the ionization balance using equation 14 in :cite:`1993A&A...279..447M`:

.. math::
        \Phi_{i,j} &= \frac{N_{i, j+1} n_e}{N_{i, j}} \\

        \Phi_{i, j} &= W \times[\delta \zeta + W ( 1 - \zeta)] \left(\frac{T_\textrm{e}}{T_\textrm{R}}\right)^{1/2}
        \Phi_{i, j}(\textrm{LTE}) \\


In the last step, we calculate the ion number densities according using the methods in :class:`LTEPlasma`

Finally we calculate the level populations (:func:`NebularPlasma.calculate_level_populations`),
by using the calculated ion species number densities:

.. math::

    N_{i, j, k}(\textrm{metastable}) &= W\frac{g_k}{Z_{i, j}}\times N_{i, j} \times e^{-\beta_\textrm{rad} E_k} \\
    N_{i, j, k}(\textrm{not metastable}) &= \frac{g_k}{Z_{i, j}}\times N_{i, j} \times e^{-\beta_\textrm{rad} E_k} \\


This concludes the calculation of the nebular plasma. In the code, the next step is calculating the :math:`\tau_\textrm{Sobolev}` using
the quantities calculated here.

Here's an example::
    l;dskfslkdfjsldkfjs

.. _tau_sobolev:

Sobolev optical depth
=====================

This function calculates the Sobolev optical depth :math:`\\tau_\\textrm{Sobolev}`



.. math::
    C_\textrm{Sobolev} = \frac{\pi e^2}{m_e c}

    \tau_\textrm{Sobolev} = C_\textrm{Sobolev}\,  \lambda\, f_{\textrm{lower}\rightarrow\textrm{upper}}\,
        t_\textrm{explosion}\, N_\textrm{lower}



.. note::
    Currently we're ignoring the term for stimulated emission:
        :math:`(1 - \frac{g_\textrm{lower}}{g_\textrm{upper}}\frac{N_\textrm{upper}}{N_\textrm{lower}})`



Macro Atom
==========

The macro atom is described in detail in :cite:`2002A&A...384..725L`. The basic principal is that when an energy packet
is absorbed that the macro atom is on a certain level. Three probabilities govern the next step the P\ :sub:`up`,
P\ :sub:`down` and P\ :sub:`down emission` being the probability for going to a higher level, a lower level and a lower
level and emitting a photon while doing this respectively (see Figure 1 in :cite:`2002A&A...384..725L` ).


The macro atom is the most complex idea to implement as a data structure. The setup is done in `~tardisatomic`, but
we will nonetheless discuss it here (as `~tardisatomic` is even less documented than this one).

For each level we look at the line list to see what transitions (upwards or downwards are possible). We create a two arrays,
the first is a long one-dimensional array containing the probabilities. Each level contains a set of probabilities, The first
part of each set contains the upwards probabilities (internal upward), the second part the downwards probabilities
(internal downward), and the last part is
  the downward and emission probability.

+--------------------------+
|Probabilities             |
+==========================+
| P\ :sub:`internal up` 1  |
+--------------------------+
| P\ :sub:`internal up` 2  |
+--------------------------+
|...                       |
+--------------------------+
| P\ :sub:`internal up` n  |
+--------------------------+
| P\ :sub:`internal down` 1|
+--------------------------+
| P\ :sub:`internal down` 2|
+--------------------------+
|...                       |
+--------------------------+
| P\ :sub:`internal down` n|
+--------------------------+
| P\ :sub:`emission down` 1|
+--------------------------+
| P\ :sub:`emission down` 2|
+--------------------------+
|...                       |
+--------------------------+
| P\ :sub:`emission down` n|
+--------------------------+

each set is stacked after the other one to make one long one dimensional `~numpy.ndarray`.

The second array is for book-keeping it has exactly the length as levels (with an example for the Si II level 15):

+--------+------------------+------------+----------------+-----------------+
|Level ID| Probability index|N\ :sub:`up`| N\ :sub:`down` | N\ :sub:`total` |
+========+==================+============+================+=================+
|14001015| ???              |17          | 5              | 17 + 5*2 = 27   |
+--------+------------------+------------+----------------+-----------------+


We now will calculate the transition probabilites, using the radiative rates in Equation 20, 21, and 22
in :cite:`2002A&A...384..725L`. Then we calculate the downward emission probability from Equation 5, the downward and
upward internal transition probabilities in :cite:`2003A&A...403..261L`.

.. math::
    p_\textrm{emission down}&= {\cal R}_{\textrm{i}\rightarrow\textrm{lower}}\,(\epsilon_\textrm{upper} - \epsilon_\textrm{lower}) / {\cal D}_{i}\\
    p_\textrm{internal down}&= {\cal R}_{\textrm{i}\rightarrow\textrm{lower}}\,\epsilon_\textrm{lower}/{\cal D}_{i}\\,
    p_\textrm{internal up}&={\cal R}_{\textrm{i}\rightarrow\textrm{upper}}\,\epsilon_{i}/{\cal D}_{i}\\,

where :math:`i` is the current level, :math:`\epsilon` is the energy of the level, and :math:`{\cal R}` is the radiative
 rates.


We ignore the probability to emit a k-packet as TARDIS only works with photon packets.
Next we calculate the radidative
rates using equation 10 in :cite:`2003A&A...403..261L`.

.. math::
    {\cal R}_{\textrm{upper}\rightarrow\textrm{lower}} &=
    A_{\textrm{upper}\rightarrow\textrm{lower}}\beta_\textrm{Sobolev}n_\textrm{upper}\\
    {\cal R}_{\textrm{lower}\rightarrow\textrm{upper}} &=
    (B_{\textrm{lower}\rightarrow\textrm{upper}}n_\textrm{lower}-
    B_{\textrm{upper}\rightarrow\textrm{lower}}n_\textrm{upper})
    \beta_\textrm{Sobolev} J_{\nu}^{b}\\,

with :math:`\beta_\textrm{Sobolev} = \frac{1}{\tau_\textrm{Sobolev}}(1-e^{-\tau_\textrm{Sobolev}})` .

using the Einstein coefficients

.. math::
    A_{\textrm{upper}\rightarrow\textrm{lower}} &= \frac{8 \nu^2 \pi^2 e^2}{m_e c}~
        \frac{g_\textrm{lower}}{g_\textrm{upper}}~f_{\textrm{lower}\rightarrow\textrm{upper}}\\
    A_{\textrm{upper}\rightarrow\textrm{lower}} &= \underbrace{\frac{4 \pi^2 e^2}{m_e c}}_{C_\textrm{Einstein}}~2 \nu^2
            \frac{g_\textrm{lower}}{g_\textrm{upper}}~f_{\textrm{lower}\rightarrow\textrm{upper}}\\
    B_{\textrm{lower}\rightarrow\textrm{upper}} &= \frac{4\pi^2 e^2}{m_e h\nu c}\,f_{\textrm{lower}\rightarrow\textrm{upper}}\\
    B_{\textrm{lower}\rightarrow\textrm{upper}} &= \underbrace{\frac{4 \pi^2 e^2}{m_e c}}_{C_\textrm{Einstein}}\frac{1}{h\nu} f_{\textrm{lower}\rightarrow\textrm{upper}}\\

    B_{\textrm{upper}\rightarrow\textrm{lower}} &= \frac{4\pi^2 e^2}{m_e h\nu c}\,f_{\textrm{lower}\rightarrow\textrm{upper}}\\
    B_{\textrm{upper}\rightarrow\textrm{lower}} &= \underbrace{\frac{4 \pi^2 e^2}{m_e c}}_{C_\textrm{Einstein}}\frac{1}{h\nu}\frac{g_\textrm{lower}}{g_\textrm{upper}}f_{\textrm{lower}\rightarrow\textrm{upper}}\\

we get

.. math::
    {\cal R}_{\textrm{upper}\rightarrow\textrm{lower}} &=
        C_\textrm{Einstein} 2 \nu^2 \frac{g_\textrm{lower}}{g_\textrm{upper}}~f_{\textrm{lower}\rightarrow\textrm{upper}}
        \beta_\textrm{Sobolev}n_\textrm{upper}\\

    {\cal R}_{\textrm{lower}\rightarrow\textrm{upper}} &=
            C_\textrm{Einstein}\frac{1}{h\nu} f_{\textrm{lower}\rightarrow\textrm{upper}}
            (n_\textrm{lower}-\frac{g_\textrm{lower}}{g_\textrm{upper}}n_\textrm{upper})
                        \beta_\textrm{Sobolev} J_{\nu}^{b}\\

This results in the transition probabilities:

.. math::
    p_\textrm{emission down}&= C_\textrm{Einstein} 2 \nu^2 \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                       \beta_\textrm{Sobolev}n_\textrm{i}\,(\epsilon_\textrm{i} - \epsilon_\textrm{lower}) / {\cal D}_{i}\\
    p_\textrm{internal down}&= C_\textrm{Einstein} 2 \nu^2 \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                       \beta_\textrm{Sobolev}n_\textrm{i}\,\epsilon_\textrm{lower}/{\cal D}_{i}\\
    p_\textrm{internal up}&=C_\textrm{Einstein}\frac{1}{h\nu} f_{\textrm{i}\rightarrow\textrm{upper}}
                                        (n_\textrm{i}-\frac{g_\textrm{i}}{g_\textrm{upper}}n_\textrm{upper})
                                                    \beta_\textrm{Sobolev} J_{\nu}^{b}\,\epsilon_{i}/{\cal D}_{i}\\,

and as we will normalise the transition probabilities numerically later,  we can get rid of :math:`C_\textrm{Einstein}`,
 :math:`\frac{1}{{\cal D}_i}` and number densities.

.. math::
    p_\textrm{emission down}&= 2 \nu^2 \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                       \beta_\textrm{Sobolev}\,(\epsilon_\textrm{i} - \epsilon_\textrm{lower})\\
    p_\textrm{internal down}&=  2 \nu^2 \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                       \beta_\textrm{Sobolev}\,\epsilon_\textrm{lower}\\
    p_\textrm{internal up}&=\frac{1}{h\nu} f_{\textrm{i}\rightarrow\textrm{upper}}
                                        \underbrace{(1-\frac{g_\textrm{i}}{g_\textrm{upper}}\frac{n_\textrm{upper}}{n_i})}
                                        _\textrm{stimulated emission}
                                                    \beta_\textrm{Sobolev} J_{\nu}^{b}\,\epsilon_{i}\\,


Like, when calculating :ref:`tau_sobolev` we assume the 'stimulated emission term' is negligible and we set it to 1.

There are two parts for each of the probabilities, one that is pre-computed by `~tardisatomic` and is in the HDF5 File,
and one that is computed during the plasma calculations:

.. math::
        p_\textrm{emission down}&= \underbrace{2 \nu^2 \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                           (\epsilon_\textrm{i} - \epsilon_\textrm{lower})}_\textrm{pre-computed}
                                           \,\beta_\textrm{Sobolev}\\
        p_\textrm{internal down} &= \underbrace{2 \nu^2 \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                           \epsilon_\textrm{lower}}_\textrm{pre-computed}\,\beta_\textrm{Sobolev}\\
        p_\textrm{internal up} &= \underbrace{\frac{1}{h\nu} f_{\textrm{i}\rightarrow\textrm{upper}}}_\textrm{pre-computed}
                                                        \beta_\textrm{Sobolev} J_{\nu}^{b}\,\epsilon_{i}\\,






.. automodapi:: tardis.plasma
