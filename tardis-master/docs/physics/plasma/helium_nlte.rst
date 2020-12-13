Helium NLTE
============

The `helium_treatment` setting in the TARDIS config. file will accept one of three options:
 * `none`: The default setting. Populate helium in the same way as the other elements.
 * `recomb-nlte`: Treats helium in NLTE using the analytical approximation outlined in an upcoming paper. 
 * `numerical-nlte`: To be implemented. Will allow the use of a separate module (not distributed with TARDIS) to perform helium NLTE calculations numerically.

Recombination He NLTE
----------------------

Paper version:

This section will summarise the equations used in the calculation of the helium state for the `recomb-NLTE` approximation in TARDIS. A full physical justification for these equations will be provided in an upcoming paper. All of the level populations are given as a function of the He II ground state population (:math:`n_{2,1,0}`), and the values are then normalised using the helium number density to give the correct level number densities.

Symbols/Indexing:
 * :math:`N_{i,j}`: Ion Number Density
 * :math:`n_{i,j,k}`: Level Number Density

  * i: atomic number
  * j: ion number
  * k: level number

.. math::
    n_{2,0,0} = 0

.. math::
    n_{2,0,k}~(k\geq1) = n_{2,1,0}\times n_{e}\times\frac{1}{W}\times\frac{g_{2,0,k}}{2g_{2,1,0}}\times\left(\frac{h^{2}}{2\pi m_{e}kT_{r}}\right)^{3/2}\times\exp{\left(\frac{\chi_{2,1}-\epsilon_{2,0,k}}{kT_{r}}}\right)\times\left(\frac{T_{r}}{T_{e}}\right)^{1/2}

(Note: An extra factor of :math:`\frac{1}{W}` is included for the metastable states of He I.)

.. math::
    n_{2,1,0} = 1

.. math::
    n_{2,1,k}~(k\geq1) = W\times\frac{g_{2,0,k}}{g_{2,1,0}}\times n_{2,1,0}\times\exp{\left(-\frac{\epsilon_{2,1,k}}{kT_{r}}\right)}

.. math::
    n_{2,2,0} = \frac{n_{2,1,0}}{n_{e}}\times[W(\delta_{2,2}\times\zeta_{2,2}+W(1-\zeta_{2,2})]\left(\frac{T_{e}}{T_{r}}\right)^{1/2}\times\frac{2g_{2,2,0}}{g_{2,1,0}}\times\left(\frac{2\pi m_{e}kT_{r}}{h^{2}}\right)^{3/2}\times\exp{\left(-\frac{\chi_{2,1}}{kT_{r}}\right)}

Code Version:

In the TARDIS plasma, some of these equations are re-written slightly to make use of existing property methods (e.g. `PhiSahaLTE`, `PhiSahaNebular`) often using the relation:

.. math::
    \frac{N_{i,j}}{Z_{i,j}} = \frac{n_{i,j,k}}{g_{i,j,k}}

Numerical He NLTE
------------------

Another `helium_treatment` option offered by TARDIS is `numerical-nlte`. The use of this plasma property requires an additional code that is the property of Stephan Hachinger (see arXiv:1201.1506) and is not distributed with TARDIS. TARDIS also requires a specific atomic datafile to use this module. This plasma option is included so that people who have access to and permission to use the necessary module may use it. Otherwise, the `recomb-NLTE` option provides a very accurate alternative approximation.