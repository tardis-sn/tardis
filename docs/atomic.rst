Atomic Data
===========
.. currentmodule:: tardis.atomic


The atomic data for tardis is stored in `hdf5 files <http://www.h5py.org/>`_. TARDIS ships with a
relatively simple atomic dataset which only contains silicon lines and levels. TARDIS also has a full atomic dataset which contains
the complete Kurucz dataset (`<http://kurucz.harvard.edu/LINELISTS/GFALL/>`_). This full dataset also contains recombination
coefficients from the ground state (:math:`\zeta-\textrm{factor}` used in :ref:`calc_zeta_label`) and data for calculating the
branching or macro atom line interaction (:ref:`macroatom`).


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
^^^^^^^^^^^
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
^^^^^^^^^^^^^^^
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
^^^^^^^^^^^
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
^^^^^^^^^^
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

