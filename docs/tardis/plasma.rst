Plasma
======
.. currentmodule:: tardis.plasma

Level Populations
-----------------

We first calculate the number densities for each atom
:math:`N_{atom} = \rho_\textrm{total} \times \textrm{Abundance fraction} / m_\textrm{atom}`

Next we calculate the partition function (Z) where i is the level number and j the ion number:
:math:`Z_{j} = \sum_{i=0}^{max levels} g_i * e^{-E_i / (k_\textrm{b} T)}`

.. automodapi:: tardis.plasma
