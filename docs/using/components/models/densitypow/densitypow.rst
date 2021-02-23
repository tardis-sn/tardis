*********************************
Using a power law density profile
*********************************

Overview
========

In this mode, the density profile (function of velocity and time since
explosion) is assumed to follow a functional form:

.. math::

    \rho (v, t_{exp}) = \rho_0 (t_{0} / t_{exp})^{3} ( v / v_{0})^{\rm exponent}

This form is defined by reference density, velocity and time parameters, and the
"exponent", each of which is set in the input file ("structure"
subsection of the "model" section, under the "density" heading; see
example below).

.. note::

   In this mode, the velocity grid has to be explicitly defined in the YAML file (see example below)


TARDIS input file example
=========================

The following example shows the relevant sections of a TARDIS input YAML file,
specifying a power law density:

.. literalinclude:: tardis_configv1_density_power_law_example.yml
    :language: yaml

