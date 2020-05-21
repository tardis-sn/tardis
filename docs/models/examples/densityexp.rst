************************************
Using an exponential density profile
************************************

Overview
========

In this mode, the density profile (function of velocity and time since
explosion) is assumed to follow a functional form:

.. math::

    \rho (v, t_{exp}) = \rho_0 (t_{0} / t_{exp})^{3} \exp( -v / v_0)

defined by reference density, velocity and time parameters. These
parameters are set in the input YAML file, specifically in the "structure"
subsection of the "model" section, under the "density" heading (see
example below).

.. note::

   In this mode, the velocity grid has to be explicitly defined in the YAML file (see example below)


TARDIS input file example
=========================

The following example shows the relevant sections of a TARDIS input YAML file,
specifying an exponential density:

.. literalinclude:: tardis_configv1_density_exponential_example.yml
    :language: yaml

