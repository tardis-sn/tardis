Atomic
======
.. currentmodule:: tardis.atomic

All atomic data files are stored in an hdf5 file

hdf5
----
All atomic data are stored in hdf5 files, in tables. The atomic data table contains information on atomic number,
symbol, mass; and the ionization data table contains atomic number, ion number, and ionization energy information.
Tables contain attributes which contains units.




Tables
------

Atomic Data

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

+------------------------+------------------------------+----------+
|Name                    | Description                  | Symbol   |
+========================+==============================+==========+
| atomic_number(z)       | Atomic Number                | 1        |
+------------------------+------------------------------+----------+
| ion number             | Ion Number                   | 1        |
+------------------------+------------------------------+----------+
| ionization energy      | Ionization Energy of atom    | eV       |
+------------------------+------------------------------+----------+




.. automodapi:: tardis.atomic