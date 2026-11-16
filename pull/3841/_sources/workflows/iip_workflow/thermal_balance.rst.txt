Type IIP thermal balance
========================

Species treatment
-----------------

H I is the only NLTE species. Its level populations are solved from the NLTE
rate equations, and the H ionization balance uses the H I continuum rates.
Elements without continuum rates, such as helium when only H I is listed under
``continuum_interaction.species``, follow the C-TARDIS treatment:

- Ionization uses the nebular approximation of Mazzali & Lucy (1993) at the
  radiation temperature and dilution factor of the current iteration and at
  each trial electron temperature. Charge conservation sets the electron
  density.
- Excitation is dilute LTE at the radiation temperature and dilution factor of
  the current iteration.

Collision strengths for levels missing from the atomic level data are ignored
with a warning, so a model with helium needs no extra atomic data.

Level-population solve
----------------------

The Type IIP workflow solves NLTE level populations for each thermal-balance
trial. When the direct root does not yield a physical population that closes
the normalized level equations to ``1e-10``, it retries from the fixed initial
population using logarithmic population ratios. These ratios keep trial
fractions positive and normalized.

If the equations remain open, the level solver retries with an empirical
``1e-4`` floor on each starting fraction and then renormalizes. This permits
excited levels depleted by the previous thermal state to recover. The solver
keeps the candidate with the smallest normalized physical equation residual,
regardless of the optimizer's reported success status.

Thermal-balance solve
---------------------

In the inner shells cooling is dominated by collisional bound-bound pairs
close to detailed balance, so the fractional-heating residual responds about
100 times less to the temperature than the electron residual does to the
density. The joint solve therefore divides each residual row by its Jacobian
diagonal, measured once per solve by forward differences at the starting
candidate. This puts both rows of every shell in candidate units without
moving the root, and lets the trust-region solver cross the curved, stiff
valley of these shells.

The thermal optimizer can evaluate an approximate physical level state during
its search. After the joint solve, the workflow rebuilds the full plasma state
from the joint candidate. That state must close the level and charge
equations, and its populations must be valid; otherwise the workflow raises an
error. A shell whose electron-density change exceeds ``2e-8`` or whose
fractional heating exceeds ``2e-7`` is reported by a warning that lists each
such shell with the ratio of its residuals to these limits, and the run
continues. The next iteration starts its thermal solve from the accepted
state, so such shells continue to converge across iterations.
