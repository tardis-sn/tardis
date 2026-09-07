The recommended end state is one configurable optical/UV Monte Carlo transport path—not separate classic, IIP, and nonhomologous solver stacks. Retain two small geometry-specific resonance kernels; consolidate everything around them. Because these APIs are internal, delete the old solver classes and modules directly without compatibility shims.

## The real physical splits

| Physical axis | Current behavior | Recommended representation |
|---|---|---|
| Kinematics/resonance | Classic and IIP use homologous expansion; nonhomologous uses a piecewise-linear velocity field and quartic resonance distance | Resolve a compatible compiled transport kernel now; reduce it to two resonance kernels after the surrounding state machines share a contract |
| Interaction processes | Classic/nonhomologous: line + electron; IIP: line + electron + bound-free/free-free | Explicit process configuration and empty continuum arrays when disabled |
| Relativity | Classic configurable; IIP forced full; nonhomologous rejects full | One explicit frame-treatment field, validated before compilation |
| Line redistribution | Scatter, downbranch, macroatom | Keep the existing enum and transition data; not a transport mode |
| Macro-atom sampling | Classic direct walk; IIP absorbing-chain sampler | A sampler choice inside shared event handling |
| Virtual packets | Classic/nonhomologous only | Optional observer/peel-off component |
| Thermal balance | Coupled to the IIP workflow | Outer iteration state update, outside packet propagation |

The current mode choice is not actually configuration driven. Standard construction hard-codes classic transport in [simulation/base.py](/home/afullard/tardis/tardis/simulation/base.py:741), while IIP and nonhomologous select their classes through dedicated workflows in [type_iip_workflow.py](/home/afullard/tardis/tardis/workflows/type_iip_workflow.py:191) and [nonhomologous_tardis_workflow.py](/home/afullard/tardis/tardis/workflows/nonhomologous_tardis_workflow.py:120).

The architectural duplication is substantial:

- Three near-copy solver classes.
- Three packet state machines.
- Separate classic/nonhomologous interaction, virtual-packet, transport-state, opacity, and radiation-field modules.
- A second IIP threaded driver whose main addition is continuum estimators. Compare the shared driver in [montecarlo_transport.py](/home/afullard/tardis/tardis/transport/montecarlo/modes/montecarlo_transport.py:238) with [iip/montecarlo_transport.py](/home/afullard/tardis/tardis/transport/montecarlo/modes/iip/montecarlo_transport.py:39).

Notably, the shared homologous tracer already accepts continuous opacity, electron-scattering probability, and a continuum-enabled flag in [homologous_rad_packet_transport.py](/home/afullard/tardis/tardis/transport/montecarlo/modes/homologous_rad_packet_transport.py:75). That is strong evidence that classic versus IIP does not merit separate transport loops.

## Recommended target

Use:

- One `MCTransportSolver`.
- One immutable resolved transport configuration.
- One `MonteCarloTransportState`.
- One Numba opacity-state type; disabled capabilities use zero-length arrays.
- One `prange` packet driver and one packet state machine.
- Two small resonance functions:
  - homologous analytic resonance, including current full-relativistic support;
  - piecewise-linear nonhomologous resonance, retaining line-direction and cursor behavior.
- Shared event functions accepting geometry or local velocity, rather than assuming `r / time_explosion`.
- One result object containing bulk, line, continuum, virtual-packet, and tracking results consistently.

The configuration should resolve these choices once before entering Numba:

```text
geometry             -> resonance kernel
continuum species    -> enabled opacity processes and estimators
line interaction     -> scatter/downbranch/macroatom
relativity setting   -> frame treatment
virtual packet count -> optional observer estimator
```

Small injected compiled functions are enough; a strategy-class hierarchy, plugin registry, or abstract solver factory would add little and complicate Numba specialization.

### Example implementation boundary

The example in `montecarlo/transport_physics.py` resolves the matching outer
driver and complete packet state machine, rather than only the line-distance
function. This is deliberately broader than the final target. Today the
nonhomologous split also owns line initialization and ordering, both line
cursors, estimator updates, virtual packets, and interaction events; injecting
only its resonance function into the classic state machine would be physically
incorrect. The resolved unit can shrink to the two resonance kernels once those
surrounding operations have a shared interface.

Relativity is resolved alongside that kernel and passed as the only
kernel-visible Boolean. Homologous time and local velocity come from the
geometry, avoiding a second kinematic source of truth. Unsupported combinations
are rejected before entering Numba.

This example preserves IIP's existing forced-full transport behavior. Packet
source selection happens earlier and still follows the input configuration, so
an IIP configuration with `enable_full_relativity: false` can pair a
nonrelativistic source with full-relativity transport. Resolving relativity
before packet-source construction would change scientific results and needs a
separate behavior change with a regression baseline.

Remove the mutable continuum global used by [montecarlo_globals.py](/home/afullard/tardis/tardis/transport/montecarlo/configuration/montecarlo_globals.py:1) and [opacity_state.py](/home/afullard/tardis/tardis/opacities/opacity_state.py:211). It makes run configuration implicit and potentially order-dependent. `OpacityStateNumba` already contains continuum fields and an absorbing-Markov field in [opacity_state_numba.py](/home/afullard/tardis/tardis/opacities/opacity_state_numba.py:13), so the separate IIP type is unnecessary.

I would also:

- Derive kinematics from the actual geometry and delete the unused `enable_nonhomologous_expansion` schema option in [montecarlo.yml](/home/afullard/tardis/tardis/io/configuration/schemas/montecarlo.yml:61).
- Derive continuum capability from configured continuum species.
- Move GPU selection entirely to the formal-integral solver; MC transport does not use it.
- Reject unsupported combinations at construction time:
  - nonhomologous + full relativity;
  - nonhomologous + continuum;
  - continuum with non-macroatom line treatment;
  - continuum virtual packets until their spawning semantics are implemented.
- Remove or explicitly reject the currently nonfunctional adiabatic-cooling and two-photon switches.

Gamma-ray transport should remain separate: it has different particles, interaction physics, deposition semantics, and energy domain. The formal integral should also remain a post-processing spectrum solver rather than joining the Monte Carlo path.

## Critical scientific boundary

Do not generalize the nonhomologous equations as part of this behavior-preserving refactor.

For a general spherical radial velocity field, the line-of-sight velocity gradient is

\[
Q(\mu)=\mu^2\frac{dv}{dr}+(1-\mu^2)\frac{v}{r},
\]

and Sobolev optical depth scales as \(1/|Q|\). Only in homologous expansion does \(dv/dr=v/r=1/t\), making this isotropic.

Current nonhomologous opacity uses only \(1/(dv/dr)\) in [tau_sobolev.py](/home/afullard/tardis/tardis/transport/montecarlo/modes/nonhomologous/tau_sobolev.py:55), while line traversal is selected from the sign of `dv/dr` in [rad_packet_transport.py](/home/afullard/tardis/tardis/transport/montecarlo/modes/nonhomologous/rad_packet_transport.py:75). TARDIS’s own documentation notes that isotropic line re-emission is special to homologous expansion in [lineinteraction.rst](/home/afullard/tardis/docs/physics_walkthrough/montecarlo/lineinteraction.rst:27).

Therefore:

- Preserve the present nonhomologous approximation exactly during consolidation.
- Keep its resonance and line-cursor logic isolated.
- Treat directional Sobolev optical depth, anisotropic escape, multiple resonances, nonhomologous continuum, and full relativity as later scientific behavior changes with separate requirements and regression baselines.

## Safest implementation order

1. Add a capability-resolution matrix and missing characterization tests. Preserve fixed-seed packet outputs and random-number draw order.

2. Consolidate the three host solvers, transport states, result contracts, tracking, threading, and configuration while retaining the existing packet functions.

3. Unify the opacity-state type and outer threaded drivers. Disabled continuum and virtual-packet paths must consume no extra random draws.

4. Merge interaction-event and virtual-packet copies by passing local velocity/geometry explicitly.

5. Merge the three packet state machines while retaining the two resonance kernels.

6. Generalize shared Sobolev and radiation-field normalization around inverse velocity gradient, verifying the homologous limit.

7. Delete the old `classic`, `iip`, and nonhomologous solver/state/event modules and migrate callers directly.

Only after equivalence is demonstrated should unsupported physical combinations be implemented.

## Validation status

Regression-backed checks were run in the `tardis` Conda environment against `/home/afullard/tardis-regression-data` at commit `631b0742fd3d2f23b8d945908c08b185ee4d71c3`:

- `test_transport.py --tardis-regression-data=…`: 25 passed.
- JIT-disabled packet dispatch selection: 9 passed.
- `test_nonhomologous.py --tardis-regression-data=…`: 26 passed.
- Classic-continuum configuration rejection: 1 passed.
- `git lfs fsck`: OK.

The regression checkout had pre-existing untracked files; they were preserved. The TARDIS worktree remains clean. No files were changed. The full classic/IIP integrations, broad suite, Ruff, and documentation build were not run because this was a read-only investigation.
