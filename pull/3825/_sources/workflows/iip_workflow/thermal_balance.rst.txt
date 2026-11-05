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

Each Monte Carlo transport before the final iteration is divided into 16
independent packet ranges. The workflow normalizes each range like the summed
estimators. Each range supplies its own continuum estimators and its own mean
intensities for the H I lines; all other lines keep the summed values. The
standard error of the electron-density and heating residuals is measured from
the scatter of the ranges at the starting candidate of the thermal solve. A
range whose ionization balance cannot be solved there is left out with a
warning. Thermal balance requires at least two packet ranges, and raises an
error otherwise.

A shell is closed when both of its residuals are within twice their standard
errors. If a 2-sigma limit exceeds ``1e-2``, the workflow warns that more
packets are needed to constrain the plasma state. The joint thermal solve
returns the starting candidate if every shell is already closed. Otherwise it
stops once every shell is closed, or once the optimizer has used its budget of
100 residual evaluations, not counting those for its finite-difference
Jacobian.

The heating row is the fractional heating ``(H - C) / C``, and its Monte
Carlo limit is measured in the same units. In the inner shells cooling is
dominated by collisional bound-bound pairs close to detailed balance, so the
heating row responds about 100 times less to the temperature than the
electron row does to the density. The joint solve therefore divides each
residual row by its Jacobian diagonal, measured once per solve by forward
differences at the starting candidate. This puts both rows of every shell in
candidate units without moving the root, and lets the trust-region solver
cross the curved, stiff valley of these shells. The limits are scaled the same
way, so the closure test is unchanged. If the ionization balance cannot be
solved at a trial point of the optimizer, that point is rejected and the
optimizer shortens its step.

If the budget runs out first, the accepted state does not close every shell.
The workflow then logs a warning that lists each shell whose electron or
heating residual is outside its 2-sigma limit, with the ratio of the residual
to that limit. The next iteration starts its thermal solve from the accepted
state, so open shells continue to converge across iterations.

Level populations must close their equations at the trial electron density
to ``1e-10`` relative to the gross rate of each equation, i.e. the sum of the
magnitudes of the rates into and out of the level. An absolute limit would
depend on the unit of time and on how fast each level is populated. Charge
conservation is solved to ``1e-10`` of the maximum electron density. The
workflow raises an error if either closure, or the validity of the
populations and Sobolev quantities, fails. Because the electron residual only
closes to its Monte Carlo limit, the trial density and the independently
charge-solved density may differ by up to that limit.
