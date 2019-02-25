***********************************************
Model with exponential density profile and uniform abundances
***********************************************

TARDIS can be used to compute a synthetic spectrum for a model with a
user-specified density and chosen set of abundances. One simple (but
flexible) class of model is when the density is parameterised as an
exponential function of velocity, as described below.


Exponential density profile
=========================

In this mode, the density profile (function of velocity and time since
explosion) is assumed to follow a functional form:

.. math::

    \rho (v, t_{exp}) = \rho_0 (t_{0} / t_{exp})^{3} \exp( -v / v_0)

defined by reference density, velocity and time parameters. These
parameters are set in the input yaml file, specifically in the "structure"
subsection of the "model" section, under the "density" heading (see
example below).


Uniform abundances
============================

For a model with an exponential density profile, a set of
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

Here is an example of an input file that sets up an exponential density
profile with a uniform set of
abundances:  

.. literalinclude:: tardis_configv1_density_exponential_test.yml
    :language: yaml

