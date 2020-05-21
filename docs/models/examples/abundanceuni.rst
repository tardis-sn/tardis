***************************
Using a uniform composition
***************************

Overview
========

The simplest possibility for specifying an ejecta composition is to use a
uniform abundance pattern in all cells specified in the density profile setup
step. These constant abundances  can be supplied directly in the input (YAML)
file. Elemental and Isotopic abundances are set in the "abundances" subsection of the "model"
section, following the "type: uniform" specifier (see example input
file below). They are specified as mass fractions, e.g.

.. code-block:: none

    Si: 0.6
    S: 0.2
    Ni56: 0.2

will set the mass fraction of silicon (Z=14) to 0.6, sulphur (Z=16) to 0.2 and Nickel (Z=26, A=56) to 0.2.

.. note::
    
    Suppose you specify Elemental and Isotopic abundance for the same element, e.g.
    
    .. code-block:: none

        Ni: 0.8
        Ni56: 0.2

    Here, Ni will refer to the stable Nickel (i.e. Z=26, A=58).

.. note::
    
    The mass fractions must sum to one. If mass fractions are supplied that do not sum to one, TARDIS will
    renormalise all the supplied abundances and print a "WARNING" message.


TARDIS input file
=================

The following example shows the relevant section of a TARDIS input YAML file
which specifies a uniform ejecta composition.

.. literalinclude:: tardis_configv1_abundance_uniform_example.yml
    :language: yaml

