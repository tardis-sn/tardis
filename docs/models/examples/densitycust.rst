.. _densitycust:

******************************
Using a custom density profile
******************************

Overview
========

TARDIS also has the capability to work with arbitrary density profiles. This is
particularly useful if the results of detailed explosion simulations should be
mapped into TARDIS. The density profile is supplied in the form of a simple
ASCII file that should look something like this:

.. literalinclude:: density.dat

In this file:

- the first line gives the reference time (see below)

- (the second line in our example is a comment)

- the remaining lines (ten in our example) give an indexed table of points that specify mass density (g / cm^3) as a function of velocity (km / s). 

TARDIS will use this table of density versus velocity to specify the density
distribution in the ejecta.  For the calculation, TARDIS will use the reference
time given in the file to scale the mass densities to whatever epoch is
requested by assuming homologous expansion:

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

Inner Boundary
==============

The first velocity-density pair in a custom density file (given by index 0) specifies the velocity of the inner boundary approximation. The density associated with this velocity is the density within the inner boundary, which does not affect TARDIS spectra. Therefore, the first density (5.4869692e-10 in the example above) can be replaced by a placeholder value. The user can choose to both specify a custom density file AND specify v_inner_boundary or v_outer_boundary in the configuration YAML file for a TARDIS run. However, the YAML-specified values must be within the velocity range specified in the custom density file, otherwise TARDIS will raise an error. When one of the YAML-specified boundary velocities falls within the velocity range specified in the custom density file, then the boundary velocity is set equal to the number in the configuration YAML file. This has the effect of splitting a cell in the custom density file into two parts, a region within the boundary and a region outside the boundary.

.. toctree::

    Custom_Density_And_Boundary_Velocities.ipynb

It is always a good idea to check the model velocities and abundances used in a TARDIS simulation after it has been successfully run.

.. warning::

   The example given here is to show the format only. It is not a
   realistic model. In any real calculation, better resolution
   (i.e. more grid points) should be used.


TARDIS input file
=================

If you create a correctly formatted density profile file (called "density.dat"
in this example), you can use it in TARDIS by putting the following lines in
the model section of the YAML file:

.. literalinclude:: tardis_configv1_density_cust_example.yml
    :language: yaml

.. note::
    The specifications for the velocities of the inner and outer boundary values can be neglected
    (in which case TARDIS will default to using the full velocity range specified in the density.txt file).
    Values for the boundary velocities that lie outside the range covered by density.txt will not be accepted.

