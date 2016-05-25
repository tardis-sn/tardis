***************************
Using a uniform composition
***************************

Overview
========

The simplest possibility for specifying an ejecta composition is to use a
uniform abundance pattern in all cells specified in the density profile setup
step. These constant abundances  can be supplied directly in the input (yaml)
file. Elemental abundances are set in the "abundances" subsection of the "model"
section, following the "type: uniform" specifier (see example input
file below). They are specified as mass fractions. E.g.

.. code-block:: none

    Si: 0.6
    S: 0.4

will set the mass fraction of silicon (Z=14) to 0.6 and sulphur (Z=16) to 0.4.

.. note::
    
    The mass fractions must sum to one. If mass fractions are supplied that do not sum to one, TARDIS will
    renormalise all the supplied abundances and print a "WARNING" message.


TARDIS input file
=================

The following example shows the relevant section of a TARDIS input yaml file
which specifies a uniform ejecta composition;

.. literalinclude:: tardis_configv1_abundance_uniform_example.yml
    :language: yaml

