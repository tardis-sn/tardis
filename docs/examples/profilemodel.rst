***********************************************
Model with custom density and abundance profile
***********************************************

TARDIS can be used to compute a synthetic spectrum for a model with user-specified density and abundance profiles -
this makes it possible to experiment with 1D profiles based on explosion models or with any empirical description of a
model with stratified abundances or density.


Arbitrary density profile
=========================

The density profile is supplied in the form of a simple ascii file that should look something like this:

.. literalinclude:: density.dat

In this file:

- the first line gives the reference time (see below)

- (the second line in our example is a comment)

- the remaining lines (ten in our example) give an indexed table of points that specify mass density (g / cm^3) as a function of velocity (km /s). 

TARDIS will use this table of density versus velocity to specify the density distribution in the ejecta.
For the calculation, TARDIS will use the reference time given in the file to scale the mass densities to whatever
epoch is requested by assuming homologous expansion:

.. math::

     \rho (t_{exp}) = \rho (t_{ref}) (t_{ref} / t_{exp})^{3}

The values in the example here define a density profile that is dropping off with

.. math::

    \rho \propto v^{-5}

.. note::

    The grid of points specified in the input file is interpreted by
    TARDIS as defining a grid in which the tabulated velocities are
    taken as the outer boundaries of grid cells and the density is
    assumed to be uniform with each cell.

.. warning::

   The example given here is to show the format only. It is not a
   realistic model. In any real calculation better resolution
   (i.e. more grid points) should be used.

Stratified abundance profile
============================

For a model with density profile supplied via a file (see above), uniform abundances can be supplied as normal.
Alternatively, a set of stratified elemental abundances can also be supplied. As with the density,
the abundances are specified via an ascii file. An ascii file that could work with the example density file given
above should be formatted like this:

.. literalinclude:: abund.dat



In this file:

- there should be the same number of rows as there were indexed points in the density profile file
- each row contains 31 numbers, the first of which is the index (i.e. matching the zone to the density profile file)
- the remaining 30 entries in each row give the set of elemental abundances for atomic number Z=1 to Z=30 (in order)

The abundances are specified as mass fractions (i.e. the sum of columns 1 to 30 in each row should be 1.0).
TARDIS does not currently include any elements heavier that Z=30.
The mass fractions specified will be adopted directly in the TARDIS calculations - so if your model is
e.g. based on an explosion simulation you may need to calculate the state of any radioactive decay chains at the correct epoch.

The example file shown here has three simple layers: 

- an innermost region (indices 0 to 2) that is composed of Si (Z=14), S (Z=16), Ar (Z=18), Ca (Z=20), Fe (Z=26), Co (Z=27) and Ni (Z=28)

- a middle region (indices 3 to 7) that is composed of O (Z=8), Mg (Z=12), Si, S, Ar and Ca

- an outer region (indices 8 and 9) that is composed of C (Z=6) and O.

.. warning::

   The example given here is to show the format only. It is not a
   realistic model. In any real calculation better resolution
   (i.e. more grid points) should be used.

.. warning::

   The calculation can be no better / more complete than the atomic
   data set. For further information on the atomic database -
   including details of how to develop your own dataset to suit your
   needs, please contact us.

TARDIS input file
=================

If you create a correctly formatted density profile file (called "density.dat") and abundance profile file (called "abund.dat"),
you can use them in TARDIS by putting the following lines in the model section of the yaml file (and remove all other lines from these sections):

.. literalinclude:: tardis_configv1_ascii_density_abund.yml
    :language: yaml

.. note::
    The specifications for the velocities of the inner and outer boundary values can be neglected
    (in which case TARDIS will default to using the full velocity range specified in the density.txt file).
    Values for the boundary velocities that lie outside the range covered by density.txt will not be accepted.

