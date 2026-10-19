Type IIP thermal balance
========================

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
