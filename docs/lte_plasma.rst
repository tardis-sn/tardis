LTE Plasma
----------

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
    >>> my_plasma = plasma.LTEPlasma({'Fe':0.5, 'Ni':0.5}, 10000, 1e-13, atom_data)
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

