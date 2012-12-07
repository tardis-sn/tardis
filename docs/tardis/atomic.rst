Atomic Data
===========
.. currentmodule:: tardis.atomic


Indices
-------

As an global index for levels we use
:math:`1000\times1000\times\textrm{Atomic Number} + 1000\times\textrm{Ion Number} + \textrm{Level Number}`
This obviously works only only for atomic data where there's less than 1000 levels. One should add a check in the
atomic data classes.


HDF5
----
All atomic data are stored in `hdf5 files <http://h5py.alfven.org/docs/intro/quick.html/>`_  which contain `Tables`_
that include mass, ionization, levels and lines data.
Tables contain attributes which contain units.




Calculating radiation field correction factors according to :cite:`1993A&A...279..447M`



Tables
------
the following data is contained in hdf5 file ``atom_data.h5``

Atomic Data
***********
dataset contains ``basic_atom_data``

+------------------------+--------------------------------+---------+
| Name                   | Description                    | Symbol  |
+========================+================================+=========+
| atomic_number          | Atomic Number (e.g. He = 2)    | z       |
+------------------------+--------------------------------+---------+
| symbol                 | Symbol (e.g. He, Fe, Ca, etc.) | None    |
+------------------------+--------------------------------+---------+
| mass                   | Average mass of atom           | u       |
+------------------------+--------------------------------+---------+


Ionization Data
***************
dataset contains ``ionization_data``

+------------------------+------------------------------+----------+
|Name                    | Description                  | Symbol   |
+========================+==============================+==========+
| atomic_number(z)       | Atomic Number                | 1        |
+------------------------+------------------------------+----------+
| ion_number             | Ion Number                   | 1        |
+------------------------+------------------------------+----------+
| ionization_energy      | Ionization Energy of atom    | eV       |
+------------------------+------------------------------+----------+


Levels Data
***********
dataset contains ``levels_data``

+------------------------+------------------------------+----------+
|Name                    | Description                  | Symbol   |
+========================+==============================+==========+
| atomic_number(z)       | Atomic Number                | 1        |
+------------------------+------------------------------+----------+
| ion_number             | Ion Number                   | 1        |
+------------------------+------------------------------+----------+
| level_number           | Level Number                 | 1        |
+------------------------+------------------------------+----------+
| energy                 | Energy of a particular level | eV       |
+------------------------+------------------------------+----------+
| g                      |                              |          |
+------------------------+------------------------------+----------+
| metastable             |                              |          |
+------------------------+------------------------------+----------+


Lines Data
**********
dataset contains ``lines_data``

+------------------------+------------------------------+----------+
|Name                    | Description                  | Symbol   |
+========================+==============================+==========+
| wavelength             | Waveslength                  |          |
+------------------------+------------------------------+----------+
| atomic_number(z)       | Atomic Number                | 1        |
+------------------------+------------------------------+----------+
| ion_number             | Ion Number                   | 1        |
+------------------------+------------------------------+----------+
| f_ul                   | Upper level probability      |          |
+------------------------+------------------------------+----------+
| f_lu                   | Lower level probability      |          |
+------------------------+------------------------------+----------+
| level_id_lower         | Upper level id               |          |
+------------------------+------------------------------+----------+
| level_id_upper         | Lower level id               |          |
+------------------------+------------------------------+----------+



Macro Atom
----------

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
|...                  |
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
 :math:`\frac{1}{{\cal D}_i}` and number densities

.. math::
    p_\textrm{emission down}&= 2 \nu^2 \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                       \beta_\textrm{Sobolev}\,(\epsilon_\textrm{i} - \epsilon_\textrm{lower})\\
    p_\textrm{internal down}&=  2 \nu^2 \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                       \beta_\textrm{Sobolev}\,\epsilon_\textrm{lower}\\
    p_\textrm{internal up}&=\frac{1}{h\nu} f_{\textrm{i}\rightarrow\textrm{upper}}
                                        (1-\frac{g_\textrm{i}}{g_\textrm{upper}}\frac{n_\textrm{upper}}{n_i})
                                                    \beta_\textrm{Sobolev} J_{\nu}^{b}\,\epsilon_{i}\\
                &= \frac{1}{h\nu} f_{\textrm{i}\rightarrow\textrm{upper}}
                                                           (1-\frac{g_\textrm{i}}{g_\textrm{upper}}
                                                            e^{-\frac{(\epsilon_\textrm{upper}-\epsilon_\textrm{i})}{k_\textrm{B} T_\textrm{rad}}})
                                                                       \beta_\textrm{Sobolev} J_{\nu}^{b}\,\epsilon_{i}\\



.. automodapi:: tardis.atomic