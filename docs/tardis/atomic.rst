Atomic
======
.. currentmodule:: tardis.atomic


hdf5
----
All atomic data are stored in `hdf5 files <http://h5py.alfven.org/docs/intro/quick.html/>`_  which contain `Tables`_
that include mass, ionization, levels and lines data.
Tables contain attributes which contain units.


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



.. automodapi:: tardis.atomic