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
    assumed to be uniform within each cell. See :ref:`model-edge-cell-quantities`.

Inner Boundary
==============

The first velocity-density pair in a custom density file (given by index 0) specifies the velocity of the inner boundary approximation. The density associated with this velocity is the density within the inner boundary, which does not affect TARDIS spectra. Therefore, the first density (5.4869692e-10 in the example above) can be replaced by a placeholder value. The user can choose to both specify a custom density file AND specify ``v_inner_boundary`` or ``v_outer_boundary`` in the configuration YAML file for a TARDIS run. The YAML-specified boundary velocity then replaces the corresponding edge of the computational domain:

- If the boundary velocity falls inside a cell of the custom density file, TARDIS truncates that cell at the boundary and discards the part of the cell outside the domain, along with any cells beyond it. The truncated cell keeps the density of the full cell, so its mass decreases.
- If the boundary velocity lies outside the velocity range of the custom density file, TARDIS issues a warning and extends the innermost or outermost cell to the boundary. The extended cell keeps its density, so its mass increases.
- ``v_inner_boundary`` must be smaller than ``v_outer_boundary``; otherwise TARDIS raises an error.

See :ref:`model-edge-cell-quantities` for details.

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
    Boundary velocities that lie outside the range covered by density.txt are accepted with a warning,
    and TARDIS extends the innermost or outermost cell to reach them (see `Inner Boundary`_).

