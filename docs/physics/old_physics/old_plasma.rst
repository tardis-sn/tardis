.. _old_plasma:

*****************
Plasma (Outdated)
*****************

.. currentmodule:: tardis.plasma_array


.. warning::
    This information (particularly code examples) is not up-to-date.


This module calculates the ionization balance and level populations in BasePlasma, give a abundance fraction, temperature
and density. After calculating the state of the plasma, these classes are able to calculate :math:`\tau_\textrm{sobolev}`
for the supernova radiative transfer. The simplest BasePlasma (:class:`BasePlasma`) only calculates the atom number densities, but serves
as a base for all BasePlasma classes. The next more complex class is `LTEPlasma` which will calculate the aforementioned quantities in
Local Thermal Equilibrium conditions (LTE). The :class:`NebularPlasma`-class inherits from `LTEPlasma` and uses a more complex
description of the BasePlasma (for details see :ref:`nebular_plasma`).

.. note::
    In this documentation we use the indices :math:`i, j, k` to mean atomic number, ion number and level number respectively.

All plasma calculations follow the same basic procedure in calculating the plasma state.
This is always accomplished with the function ``update_radiationfield``. This block diagram shows the basic procedure



Base Plasma
-----------

`BasePlasma` serves as the base class for all plasmas and can just calculate the atom number densities for a given input of
abundance fraction.

.. math::
    N_{atom} = \rho_\textrm{total} \times \textrm{Abundance fraction} / m_\textrm{atom}

In the next step the line and level tables are purged of entries that are not represented in the
abundance fractions are saved in `BasePlasma.levels` and `BasePlasma.lines`. Finally, the function `BasePlasma.update_t_rad` is called
at the end of initialization to update the plasma conditions to a new :math:`T_\textrm{radiation field}` (with the give t_rad).
This function is the same in the other plasma classes and does the main part of the calculation. In the case of `BasePlasma` this is only
setting `BasePlasma.beta_rad` to :math:`\frac{1}{k_\textrm{B}T_\textrm{rad}}`.


Here's an example how to instantiate a simple base plasma::


    >>> from tardis import plasma
    >>> from tardis.io import atomic
    >>> atom_data = atomic.AtomData.from_hdf()
    >>> my_plasma = plasma.BasePlasma({'Fe':0.5, 'Ni':0.5}, 10000, 1e-13, atom_data)
    >>> print my_plasma.abundances
    atomic_number abundance_fraction number_density
    ------------- ------------------ --------------
               28                0.5  513016973.936
               26                0.5  539183641.472

Plasma Types
------------
.. toctree::
    :maxdepth: 0

    ../plasma/lte_plasma.rst
    ../plasma/nebular_plasma.rst


.. _tau_sobolev:


Sobolev optical depth
---------------------

This function calculates the Sobolev optical depth :math:`\tau_\textrm{Sobolev}`



.. math::
    C_\textrm{Sobolev} = \frac{\pi e^2}{m_e c}

    \tau_\textrm{Sobolev} = C_\textrm{Sobolev}\,  \lambda\, f_{\textrm{lower}\rightarrow\textrm{upper}}\,
        t_\textrm{explosion}\, N_\textrm{lower}
        (1 - \frac{g_\textrm{lower}}{g_\textrm{upper}}\frac{N_\textrm{upper}}{N_\textrm{lower}})




.. toctree::
    :maxdepth: 1

    ../plasma/macroatom.rst
    ../plasma/nlte.rst


.. .. automodapi:: tardis.plasma_array

