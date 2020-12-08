.. _atomic-data-descritpion:

******************************
Atomic Data Format Description
******************************
.. currentmodule:: tardis.io.atomic


The atomic data for tardis is stored in `HDF5 files <http://www.h5py.org/>`_. TARDIS ships with a
relatively simple atomic dataset that only contains silicon lines and levels. TARDIS also has a full atomic dataset which contains
the complete Kurucz dataset (`<http://kurucz.harvard.edu/LINELISTS/GFALL/>`_). This full dataset also contains recombination
coefficients from the ground state (:math:`\zeta-\textrm{factor}` used in :ref:`calc_zeta_label`) and data for calculating the
branching or macro atom line interaction (:ref:`macroatom`).



HDF5 Dataset
------------

As mentioned previously, all atomic data is stored in `HDF5 files <http://www.h5py.org/>`_  that contain tables
that include mass, ionization, levels and lines data. The atom data that ships with TARDIS is located in data/atom.

The dataset ``basic_atom_set`` contains the Atomic Number, Symbol of the elements and average mass of the elements.

Basic Atomic Data
^^^^^^^^^^^^^^^^^

+------------------------+--------------------------------+---------+
| Name                   | Description                    | Unit    |
+========================+================================+=========+
| atomic_number          | Atomic Number (e.g. He = 2)    | z       |
+------------------------+--------------------------------+---------+
| symbol                 | Symbol (e.g. He, Fe, Ca, etc.) | None    |
+------------------------+--------------------------------+---------+
| mass                   | Average mass of atom           | u       |
+------------------------+--------------------------------+---------+


The ionization data is stored in ``ionization_data``.

Ionization Data
^^^^^^^^^^^^^^^

+------------------------+------------------------------+----------+
|Name                    | Description                  | Unit     |
+========================+==============================+==========+
| atomic_number(z)       | Atomic Number                | 1        |
+------------------------+------------------------------+----------+
| ion_number             | Ion Number                   | 1        |
+------------------------+------------------------------+----------+
| ionization_energy      | Ionization Energy of atom    | eV       |
+------------------------+------------------------------+----------+


.. note:: In TARDIS, Ion 0 is neutral. 
The levels data is stored in ``levels_data``.

Levels Data
^^^^^^^^^^^


+------------------------+------------------------------+----------+
|Name                    | Description                  | Unit     |
+========================+==============================+==========+
| atomic_number(z)       | Atomic Number                | 1        |
+------------------------+------------------------------+----------+
| ion_number             | Ion Number                   | 1        |
+------------------------+------------------------------+----------+
| level_number           | Level Number                 | 1        |
+------------------------+------------------------------+----------+
| energy                 | Energy of a particular level | eV       |
+------------------------+------------------------------+----------+
| g                      |                              | 1        |
+------------------------+------------------------------+----------+
| metastable             |                              | bool     |
+------------------------+------------------------------+----------+

All lines are stored in ``lines_data``.

Lines Data
^^^^^^^^^^

+------------------------+------------------------------+----------+
|Name                    | Description                  | Unit     |
+========================+==============================+==========+
| wavelength             | Waveslength                  | angstrom |
+------------------------+------------------------------+----------+
| atomic_number(z)       | Atomic Number                | 1        |
+------------------------+------------------------------+----------+
| ion_number             | Ion Number                   | 1        |
+------------------------+------------------------------+----------+
| f_ul                   | Upper level probability      | 1        |
+------------------------+------------------------------+----------+
| f_lu                   | Lower level probability      | 1        |
+------------------------+------------------------------+----------+
| level_id_lower         | Upper level id               | 1        |
+------------------------+------------------------------+----------+
| level_id_upper         | Lower level id               | 1        |
+------------------------+------------------------------+----------+

The next three datasets are only contained in the full dataset available upon request from the authors.

The factor correcting for photo-ionization from excited levels (needed in :ref:`calc_zeta_label`) is stored in the dataset ``zeta_data``.
The data is stored in a special way as one large :py:class:`numpy.ndarray` where the first two columns are Atomic Number and Ion
Number. All further columns are the :math:`\zeta-\textrm{factors}` for different temperatures. The temperatures are stored
in the attribute ``t_rads``.

+------------------------+------------------------------+----------+
|Name                    | Description                  | Unit     |
+========================+==============================+==========+
| atomic_number(z)       | Atomic Number                | 1        |
+------------------------+------------------------------+----------+
| ion_number             | Ion Number                   | 1        |
+------------------------+------------------------------+----------+
| T_XXXX                 | Temperature for column       | K        |
+------------------------+------------------------------+----------+
| ...                    |  ...                         | ...      |
+------------------------+------------------------------+----------+
| T_XXXX                 | Temperature for column       | K        |
+------------------------+------------------------------+----------+



There are two datasets for using the macro atom and branching line interactions. The ``macro_atom_data`` and ``macro_atom_references``:


The ``macro_atom_data`` contains blocks of transition probabilities, several indices and flags. The Transition Type flag
has three states:

* -1 for downwards emitting
* 0 for downwards internal
* 1 for upwards internal (for more explanations, please
refer to :ref:`macroatom`).

Macro Atom Data
^^^^^^^^^^^^^^^
+-------------------------+------------------------------+----------+
|Name                     | Description                  | Unit     |
+=========================+==============================+==========+
| atomic_number(z)        | Atomic Number                | 1        |
+-------------------------+------------------------------+----------+
| ion_number              | Ion Number                   | 1        |
+-------------------------+------------------------------+----------+
| source_level_number     | Source Level Number          | 1        |
+-------------------------+------------------------------+----------+
| destination_level_number| Destination Level Number     | 1        |
+-------------------------+------------------------------+----------+
| transition_type         | Transition Type              | 1        |
+-------------------------+------------------------------+----------+
| transition_probability  | Transition Probability       | 1        |
+-------------------------+------------------------------+----------+
| transition_line_id      | Transition Line ID           | 1        |
+-------------------------+------------------------------+----------+

Here's the structure of the probability block. The atomic number, ion number and source level number are the same
within each block, the destination level number the transition type and transition probability are changing.
The transition probabilities are only part of the final probability and will be changed during the calculation.
For details on the macro atom, please refer to :ref:`macroatom`.

+-------------+------------+-------------------+------------------------+-----------------+--------------------------+--------------------+
|Atomic Number|Ion Number  |Source Level Number|Destination Level Number| Transition Type |Transition probabilities  |Transition Line ID  |
+=============+============+===================+========================+=================+==========================+====================+
| Z\ :sub:`1` | I\ :sub:`1`| i\ :sub:`1`       | j\ :sub:`1`            |        -1       | P\ :sub:`emission down` 1|       k\ :sub:`1`  |
+-------------+------------+-------------------+------------------------+-----------------+--------------------------+--------------------+
| Z\ :sub:`1` | I\ :sub:`1`| i\ :sub:`1`       | j\ :sub:`2`            |        -1       | P\ :sub:`emission down` 2|       k\ :sub:`2`  |
+-------------+------------+-------------------+------------------------+-----------------+--------------------------+--------------------+
| ...         | ...        | ...               | ...                    |        ...      | ...                      |          ...       |
+-------------+------------+-------------------+------------------------+-----------------+--------------------------+--------------------+
| Z\ :sub:`1` | I\ :sub:`1`| i\ :sub:`1`       | j\ :sub:`n`            |        -1       | P\ :sub:`emission down` n|       k\ :sub:`n`  |
+-------------+------------+-------------------+------------------------+-----------------+--------------------------+--------------------+
| Z\ :sub:`1` | I\ :sub:`1`| i\ :sub:`1`       | j\ :sub:`1`            |        0        | P\ :sub:`internal down` 1|       k\ :sub:`1`  |
+-------------+------------+-------------------+------------------------+-----------------+--------------------------+--------------------+
| Z\ :sub:`1` | I\ :sub:`1`| i\ :sub:`1`       | j\ :sub:`2`            |        0        | P\ :sub:`internal down` 2|       k\ :sub:`2`  |
+-------------+------------+-------------------+------------------------+-----------------+--------------------------+--------------------+
| ...         | ...        | ...               | ...                    |        ...      | ...                      |          ...       |
+-------------+------------+-------------------+------------------------+-----------------+--------------------------+--------------------+
| Z\ :sub:`1` | I\ :sub:`1`| i\ :sub:`1`       | j\ :sub:`n`            |        0        | P\ :sub:`internal down` n|       k\ :sub:`n`  |
+-------------+------------+-------------------+------------------------+-----------------+--------------------------+--------------------+
| Z\ :sub:`1` | I\ :sub:`1`| i\ :sub:`1`       | j\ :sub:`1`            |        1        | P\ :sub:`internal up`   1|       k\ :sub:`1`  |
+-------------+------------+-------------------+------------------------+-----------------+--------------------------+--------------------+
| Z\ :sub:`1` | I\ :sub:`1`| i\ :sub:`1`       | j\ :sub:`2`            |        1        | P\ :sub:`internal up`   2|       k\ :sub:`2`  |
+-------------+------------+-------------------+------------------------+-----------------+--------------------------+--------------------+
| ...         | ...        | ...               | ...                    |        ...      | ...                      |          ...       |
+-------------+------------+-------------------+------------------------+-----------------+--------------------------+--------------------+
| Z\ :sub:`1` | I\ :sub:`1`| i\ :sub:`1`       | j\ :sub:`n`            |        1        | P\ :sub:`internal up`   n|       k\ :sub:`n`  |
+-------------+------------+-------------------+------------------------+-----------------+--------------------------+--------------------+

The ``macro_references`` dataset contains the numbers for each block:

Macro Atom References
^^^^^^^^^^^^^^^^^^^^^

+-------------------------+------------------------------+----------+
|Name                     | Description                  | Unit     |
+=========================+==============================+==========+
| atomic_number(z)        | Atomic Number                | 1        |
+-------------------------+------------------------------+----------+
| ion_number              | Ion Number                   | 1        |
+-------------------------+------------------------------+----------+
| source_level_number     | Source Level Number          | 1        |
+-------------------------+------------------------------+----------+
| count_down              | Number of down transitions   | 1        |
+-------------------------+------------------------------+----------+
| count_up                | Number of up transitions     | 1        |
+-------------------------+------------------------------+----------+
| count_total             | Total number of transitions  | 1        |
+-------------------------+------------------------------+----------+

The Atom Data Class
-------------------

Atom Data is stored inside TARDIS in the :class:`AtomData`-class. The class method :func:`AtomData.from_hdf` will
instantiate a new `AtomData`-class from an HDF5 file. If none is given it will automatically
take the default HDF5-dataset shipped with TARDIS. A second function :func:`AtomData.prepare_atom_data`
will cut the levels and lines data to only the required atoms and ions. In addition, it will create the intricate system
of references needed by macro atom or branching line interactions.


Indexing fun
------------

The main problem with the atomic data is indexing. Most of these references require multiple numbers, e.g. atomic number,
ion number and level number. The `pandas <https://pandas.pydata.org/>`_ framework provides the ideal functions to accomplish this. In TARDIS, we extensively
use :py:class:`pandas.MultiIndex`, :py:class:`pandas.Series` and :py:class:`pandas.DataFrame`.

TO BE BETTER DOCUMENTED ...





.. automodapi:: tardis.io.atomic
