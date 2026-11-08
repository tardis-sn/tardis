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

The thermal optimizer can evaluate an approximate physical level state during
its search. If its joint-shell candidate fails to converge, the workflow
retries each shell with an open electron or fractional-heating residual,
starting from that candidate. A stalled local solve is retried from the same
shell's pre-optimization state. The workflow then rebuilds the full plasma
state. The accepted state must pass the final level, electron, charge, and
heating convergence checks; otherwise the workflow raises an error.

The local objective retains the joint optimizer's fixed level-population guess.
The final rebuild uses the resulting population as its guess, so a change to a
different nonlinear level branch can still fail final convergence.
