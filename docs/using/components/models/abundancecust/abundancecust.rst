.. _abundancecust:

*************************************
Using a custom stratified composition
*************************************

ASCII Format
=============

To use a stratified ejecta composition in TARDIS, the elemental abundances may
be specified on a per-cell basis via an external ASCII file (similar to setting
up a :ref:`custom density <densitycust>` profile). An ASCII file that could
work on a mesh with ten cells should be formatted like this:

.. literalinclude:: abund.dat

In this file:

- there should be the same number of rows as there are cells specified in the velocity/density structure part of the TARDIS setup
- each row contains 31 numbers, the first of which is the index (i.e. matching the zone to the density profile file)
- the remaining 30 entries in each row give the set of elemental abundances for atomic number Z=1 to Z=30 (in order)

The abundances are specified as mass fractions (i.e. the sum of columns 1 to 30
in each row should be 1.0). TARDIS does not currently include any elements
heavier that Z=30. The mass fractions specified will be adopted directly in
the TARDIS calculations --- so if your model is, for example, based on an explosion
simulation, you may need to calculate the state of any radioactive decay chains
at the correct epoch.

The example file shown here has three simple layers: 

- an innermost region (indices 0 to 2) that is composed of Si (Z=14), S (Z=16), Ar (Z=18), Ca (Z=20), Fe (Z=26), Co (Z=27) and Ni (Z=28)

- a middle region (indices 3 to 7) that is composed of O (Z=8), Mg (Z=12), Si, S, Ar and Ca

- an outer region (indices 8 and 9) that is composed of C (Z=6) and O.

.. warning::

   The example given here is to show the format only. It is not a
   realistic model. In any real calculation, better resolution
   (i.e. more grid points) should be used.

.. warning::

   The calculation can be no better / more complete than the atomic
   data set. For further information on the atomic database ---
   including details of how to develop your own dataset to suit your
   needs, please contact us.

TARDIS ascii input file
=======================

If you create a correctly formatted abundance profile file (called "abund.dat"
in this example), you can use it in TARDIS by putting the following lines in
the model section of the YAML file:

.. literalinclude:: tardis_configv1_abundance_cust_example.yml
    :language: yaml


CSV Format
==========

In this format, both elemental and isotopic abundances may
be specified on a per-cell basis via an external CSV file. A CSV file that could
work on a mesh with ten cells should be formatted like this:

.. literalinclude:: tardis_model_abund.csv

.. note::
    
    The file should always contain the cell index as a running index in the
    first column.

.. danger::

    The header line for the isotopic abundance structure can under no
    circumstances start with a '#'.

In this file:

- Header row contains element and isotopes symbol 
- the remaining entries in each row give the set of elemental and isotopic abundances.
- the first column contains a running index

The abundances are specified as mass fractions (i.e. the sum of columns
in each row should be 1.0). The mass fractions specified will be adopted directly in
the TARDIS calculations.

The example file shown here has three simple layers: 

- an innermost region that is composed of Si and two Nickel Isotopes Ni56 and Ni58

- a middle region that is composed of O and Mg

- an outer region that is composed of C and O

.. note::
    
    Suppose you specify Elemental and Isotopic abundances for the same element. For ex-
    :code:`Ni` and :code:`Ni56`. 
    Here, Ni will refer to the stable Nickel (i.e. Z=26, A=58).


.. note::
  
    As with the custom density file, the first row will be ignored. It is
    supposed to give the composition below the photosphere. Thus, the first row
    (after the header) can be filled with dummy values.

.. warning::

   The example given here is to show the format only. It is not a
   realistic model. In any real calculation better resolution
   (i.e. more grid points) should be used.

TARDIS CSV input file
=====================

If you create a correctly formatted abundance profile file (called "tardis_model_abund.csv"
in this example), you can use it in TARDIS by putting the following lines in
the model section of the YAML file:

.. literalinclude:: tardis_configv1_isotope_abundance_cust_example.yml
    :language: yaml   

Convert ASCII abundance file format to CSV format
=================================================

If you want to convert an ASCII abundance file (say "abund.dat") to CSV format, you can use 
:code:`convert_abundances_format` function for it. Here is an example to demonstrate this:

.. code:: python

    from tardis.util import convert_abundances_format
    df = convert_abundances_format('abund.dat')
    df.to_csv('converted_abund.csv', index=False, sep=' ')  
