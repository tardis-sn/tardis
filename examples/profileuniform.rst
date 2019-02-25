***********************************************
Model with custom density profile and uniform abundances
***********************************************

TARDIS can be used to compute a synthetic spectrum for a model with a
user-specified density and chosen set of abundances -
this makes it possible to experiment with 1D density profiles and
simple specifications of the composition


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

Uniform abundances
============================

For a model with density profile supplied via a file (see above),
uniform abundances can be supplied directly in the input (yaml)
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

If you create a correctly formatted density profile file (called
"density.dat"), here is an example of how it can be combined with a uniform set of
abundances:  

.. literalinclude:: tardis_configv1_ascii_density_uniabund.yml
    :language: yaml

.. note::

    The specifications for the velocities of the inner and outer boundary values can be neglected
    (in which case TARDIS will default to using the full velocity range specified in the density.txt file).
    Values for the boundary velocities that lie outside the range covered by density.txt will not be accepted.

