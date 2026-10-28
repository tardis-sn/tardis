# Nonhomologous line search: traversal state implementation plan

## Science requirement and scope

A real packet must encounter a given line at most once while its projected
velocity changes monotonically along a free-flight trajectory. A line may be
encountered again after a genuine reversal of the projected velocity gradient.
Count a resonance as traversed once its Sobolev optical depth and line
estimators have been processed, whether or not it causes a scattering event.

This is a plan revision in response to
[the static audit](nonhomologous_line_search_plan_audit.md); no transport code
or tests have been changed yet. It covers the real-packet line search in
`tardis/transport/montecarlo/modes/nonhomologous/rad_packet_transport.py` and
the line-range helpers in `tardis/transport/geometry/calculate_distances.py`.
The packet loop in `modes/nonhomologous/packet_propagation.py` must carry search
state between calls to `trace_packet`. The separate nonhomologous virtual-packet
tracer uses the sign of radial `dv/dr` and does not split at projected-gradient
zeros. Assess that path separately. Until then, this fix establishes
real-packet line-search behavior; it does not establish corrected spectra.

Assumptions: shell velocity is piecewise linear in radius; the transport line
list is sorted by descending rest frequency; nonhomologous real-packet
transport uses its existing first-order Doppler treatment. Preserve the
current exception for a resonance at exactly zero projected gradient, where
the Sobolev optical depth is singular. No issue identifier was supplied.

## Reproduced defect

For one shell with radii `[1e14, 2e14]` cm and edge velocities
`[1e8, 1.8e8]` cm/s, take a packet at `r = 1.25e14` cm, `mu = 0.1`, and a
line at `1e15` Hz. Set the lab-frame packet frequency to
`nu_line / (1 - v(r) * mu / c)`. The field is genuinely nonhomologous:
`v(r) = (8e-7 s^-1) r + 2e7 cm/s`. Its projected gradient stays positive
along this ray, so there is no physical opportunity to revisit the line.

As a direct diagnostic, set `next_line_id = 1` and `prev_line_id = 0` to
indicate that this line was crossed. Current `trace_packet` still selects it
about **1.455 cm** ahead and adds a line-estimator contribution. This packet-ID
setup is not the primary red test: those IDs will not be the new traversal
state.

The primary red case uses two actual traces. In the same shell, launch from
`r = 1.25e14` cm and `mu = 0.1`. Construct the lab frequency so the first
resonance of the `1e15` Hz line lies `1e12` cm along the ray. Give it high line
optical depth and negligible electron opacity. After the first `LINE` event,
move the packet to that resonance, set the coherent outgoing direction to
`mu = 1`, and emit at the same line. On the current implementation the second
trace returns another `LINE` event about **1.891 cm** ahead, instead of the
outer boundary about `7.49e13` cm away, and increments the estimator again.
Both traces have positive projected gradient. The shell-width-scaled distance
cutoff is **1 cm**. The exact short distance may vary with floating-point and
compiler settings; the event sequence is the test contract.

The mechanism has two parts:

1. `get_line_id_range_nonhomologous` includes frequencies near both interval
   endpoints. Each call to `trace_packet` rebuilds this range from the packet's
   comoving frequency and ignores prior traversal state.
2. `calculate_distance_line_nonhomologous` subtracts large ray coordinates to
   recover a small distance. Rounding can put a previously crossed resonance
   slightly ahead of the packet. Its fixed distance cutoff cannot reliably
   distinguish that artifact from a distinct nearby resonance.

The existing `next_line_id` and `prev_line_id` cannot serve as the traversal
record: frequency-based updates overwrite them, and the interaction handlers
use them to identify absorbed or emitted lines.

## Mathematical basis

For `v(r) = m r + q`, let `x = r_0 mu_0 + s` be signed distance along the ray
and `b^2 = r_0^2 (1 - mu_0^2)` its squared impact parameter. Then

```text
r(x) = sqrt(x^2 + b^2)
w(x) = v(r) x/r = m x + q x/sqrt(x^2 + b^2)
G(x) = dw/ds = m + q b^2/(x^2 + b^2)^(3/2)
d nu_comov/ds = -(nu_lab/c) G(x)
```

Genuine, resolved zeros of `G` partition a flight into intervals with one
line-search direction. Where `G > 0`, comoving frequency decreases and line IDs
increase; where `G < 0`, frequency increases and IDs decrease. A zero that
only touches `G = 0` without changing its sign does not create a new interval.
The line IDs visited within a monotonic interval have a contiguous directional
frontier, so one cursor suffices for that interval.

The current zero finder drops roots within `CLOSE_LINE_THRESHOLD * shell_width`
of the ray start or shell edge and suppresses a near-tangent pair using the
same *radial* tolerance. That radial test is not a bound on separation along
the ray: near tangency, `x_zero = sqrt(r_zero^2 - b^2)` can be much larger than
`r_zero - b`. Before claiming complete interval coverage, revise or test this
filter. Retain sign-changing zeros at distinct representable forward ray
coordinates, even when the old radial or endpoint-distance cutoff would omit
them. A zero that rounds onto the ray start, shell edge, or its paired zero
cannot define a separate floating-point interval. Record that coordinate
resolution limit explicitly. Also bound the comoving-frequency excursion
between each endpoint and extremum by propagating floating-point error through
`w(x)` and `nu_lab (1 - w/c)`; include cancellation in the two terms of `w`,
not only spacings of the final frequency. If that excursion is within the
bound, label line re-encounters there spectrally unresolved rather than
claiming correctness. Compare forward ray distances with the existing solver
cutoff separately; a radial distance must not stand in for that check. Test
near-edge and near-tangent cases on both sides of these bounds. Keep the
ordinary two-root physics test well away from `G = 0`.

## Proposed search state and algorithm

Carry two Numba-compatible integers owned by the nonhomologous packet loop:
`next_untraversed_line_id` and `line_search_direction`. A direction of `0`
means the cursor has no active frontier and must be initialized from the next
interval's frequency range. A direction of `+1` or `-1` means the cursor is
the first unprocessed line ID in that direction; `line_count` and `-1` are
valid exhausted cursors. This invariant needs no separate validity flag.
Initialize direction to `0` for a newly launched packet. Add the two scalar
state values to `trace_packet`'s arguments and return the updated values with
every `LINE`, `ESCATTERING`, or `BOUNDARY` result. `packet_propagation` retains
them between loop iterations. This avoids a new mutable Numba class and keeps
the public packet fields available for their current interaction semantics.
Update direct tracer callers in tests for the new return contract. Keep this
state separate from packet IDs used by line-interaction code.

For each analytic interval:

1. Determine the traversal direction from the analytic `G` at a nonzero
   interior point. Change `get_line_id_range_nonhomologous` to accept that
   `line_id_step` explicitly. The helper computes unordered frequency bounds
   from the two endpoint frequencies, applies its existing frequency padding,
   finds the candidate ID interval, and returns it in the supplied direction.
   It must not infer direction from `comov_nu_end < comov_nu_start`.
2. If the active direction matches `line_search_direction`, intersect the
   directional ID range with the cursor. For redshift, start at least at
   `next_untraversed_line_id`; for blueshift, start at most there. If direction
   is `0` or has reversed, start from the frequency-derived range.
3. Solve distances only for eligible lines. Advance the cursor to
   `line_id + direction` immediately after processing a real resonance, before
   any return for a line event. Do not mark a candidate if electron scattering
   or a boundary occurs first.
4. When an analytic zero reverses `G`, clear the old frontier and derive a new
   range for the next interval. A line at that zero retains the existing
   singular Sobolev-depth exception; it is not assigned to an ordinary
   interval endpoint.

For equal or nearly equal rounded endpoint frequencies with known nonzero
`G`, the helper still returns the candidate IDs in the supplied direction.
Add a focused helper test for both signs and a physical short-interval test.
This is needed for correctness near a turning point.

Apply the same sign comparison when a new trace begins after an event or shell
crossing. A shell edge can reverse `G` discontinuously even without an interior
analytic zero. Preserve the cursor across a shell edge only when its stored
outgoing and incoming edge velocities are exactly equal in the represented model,
the packet direction and lab frequency are continuous, and the one-sided
gradient signs agree. `Radial1DGeometry` permits approximate velocity equality
and `NumbaRadial1DGeometry` accepts edge arrays directly, so a permitted
velocity jump must clear the cursor: comoving frequency then jumps and can
skip a line without a continuous Sobolev resonance. Do not add broad input
validation here or claim that the reset solves transport through a velocity
discontinuity. Use stored edge values for this comparison rather than two
independent boundary interpolations, which may round differently. An outgoing
packet created at a line resonance must start
beyond that emitted line, even if its new direction has a different gradient
sign, to avoid immediately replaying the emission point.

| Event | Traversal-state action |
| --- | --- |
| Line crossed without scattering | Advance past the processed line. |
| Coherent line scattering | Keep the frontier if the outgoing gradient sign agrees; otherwise begin a new directional interval. Exclude the emitted line at its starting point. |
| Thomson scattering | Preserve the frontier for unchanged gradient sign; reset it if the new ray reverses that sign. |
| Emission into a different line through macro-atom or downbranch processing | Begin a new frequency history and seed the emitted line as traversed. The earlier directional range is no longer contiguous relative to the new frequency. |
| Shell boundary | Preserve only for exact represented velocity continuity and matching one-sided gradient signs. Otherwise clear the frontier; a velocity jump is outside the continuous-resonance claim. |
| Interior analytic gradient reversal | Clear the frontier before searching the next interval. |

After a line event, compare the absorbed line ID recorded by `trace_packet`
with the emitted line ID set by `line_emission`; use IDs rather than a
floating-point frequency comparison to detect a new emission channel. Update
the outgoing cursor after `line_scatter_event` has chosen the new direction
and frequency, before calling `trace_packet` again. After Thomson scattering,
compare the incoming and outgoing projected-gradient signs while retaining
the cursor for an elastic, same-sign event.

### Endpoint ownership

An exact resonance at a shell edge belongs to the **upstream shell along the
packet path**. Use that shell's Sobolev coefficient, projected gradient, and
line estimators for both outward and inward packets. This takes the one-sided
Sobolev limit from the region the packet is leaving. Process the line before
the boundary event, then carry the advanced cursor into the adjacent shell
when velocity and gradient direction are continuous. Define the edge case by
the unsquared resonance equation evaluated at the boundary; admit a solver
root displaced by demonstrable floating-point error only when its residual
is compatible with that equation. Clamp an admitted edge event to the boundary
distance. A broad frequency padding alone must not turn a merely nearby line
into an edge event. The downstream shell treats the start point as excluded.
Test one continuous two-shell field with different shell opacities and assert
one upstream estimator contribution and no downstream duplicate. Parameterize
inward and outward travel if practical. If the upstream projected gradient is
exactly zero at an edge resonance, retain the singular-depth exception.

Keep the existing distance tolerance for numerical root filtering. The cursor
provides the identity-based guard against replay; increasing a spatial
epsilon could suppress a physically distinct close line. The current solver
already rejects roots inside its `1e-14 * shell_width` spatial cutoff, so this
change does not promise resolution of arbitrarily close physical lines.

## Test-first implementation sequence

1. Write `test_nonhomologous_coherent_line_does_not_replay` against the
   existing packet-loop interface, with the real tracer and real packet move.
   Use the two-trace scenario above, a controlled `mu = 1` coherent outgoing
   direction, one high-opacity line, and negligible electron opacity. Fix the
   angle draw at its sampling boundary or use an equivalent deterministic
   fixture; do not mock `trace_packet`. Confirm that the current code fails
   because it records a second `LINE` event or estimator contribution after
   the first move and emission. The implemented path must instead reach the
   analytically calculated outer boundary with exactly one line contribution.
   Bound the test's event recording so the pre-change replay cannot spin
   indefinitely; record real tracer results rather than replacing them.
   The existing `test_nonhomologous_packet_propagation_dispatch_numba_disabled`
   mocks transport dispatch and cannot establish this claim. A direct tracer
   test that injects explicit cursor state may supplement the packet-loop test
   after the cursor contract exists; an API mismatch is not a valid red phase.
2. Add `test_nonhomologous_line_reappears_after_gradient_reversal` with the
   following independent ray construction. Set `b = 1e14` cm, launch at
   `r = b`, `mu = 0`, and use a shell from `0.9e14` to `1.7e14` cm. Let
   `x_1 = 0.25 b`, `x_2 = 1.25 b`, and `q = 2e7` cm/s. Choose
   `m = q (x_2/r_2 - x_1/r_1)/(x_1 - x_2)`, where
   `r_i = sqrt(x_i^2 + b^2)`. This gives `m = -1.0766663688133946e-7 s^-1`
   and `w(x_1) = w(x_2) = 2.159046578693174e6` cm/s. For
   `nu_lab = 1e15` Hz, set `nu_line = nu_lab (1 - w(x_1)/c)`, about
   `9.999279819581488e14` Hz. The unsquared equation independently gives
   resonance distances `s_1 = 2.5e13` and `s_2 = 1.25e14` cm. The analytic
   gradient zero is at `s = 7.149272474576786e13` cm, strictly between
   them; `G(x_1) = 7.494842196954653e-8 s^-1` and
   `G(x_2) = -5.8909950242945375e-8 s^-1`. Disable line scattering and use
   negligible electron opacity so both crossings are independent of
   stochastic event selection. Assert both solver positions and the trace estimator sum from
   these independently calculated gradients: for packet energy `E`,
   `energy_deposition_line_rate = 2 E nu_line/nu_lab` and
   `mean_intensity_blueward = E (nu_line/nu_lab)/nu_lab
   * (1/abs(G_1) + 1/abs(G_2))`. For `E = 1`, these are about
   `1.9998559639162976` and `3.031538667815347e-8`, respectively.
3. Test the range helper with equal and nearly equal endpoint frequencies and
   an independently known nonzero gradient, in both directions. Assert the
   returned ID order and that a physical line within the frequency bounds is
   not skipped. Keep this helper-contract test separate from the reversal
   physics test. Do not count a missing new helper argument as a red-phase
   failure; test the behavior after its signature exists, and falsify it by
   temporarily restoring endpoint-derived direction.
4. Add a continuous two-shell edge test with a line exactly at the shared
   boundary, high upstream line opacity, zero downstream line opacity, and
   nonzero upstream `G`. With scattering enabled, assert a `LINE` event and
   upstream estimator contribution at the edge. With scattering disabled,
   assert a subsequent boundary event and no downstream duplicate. Add an
   inward counterpart if practical.
   Test same-sign preservation using exactly shared velocity edges, and test
   sign-change reset separately. For an allowed velocity jump, assert only
   that the cursor clears; do not claim a continuous-resonance result.
5. Add narrow event tests with an actual line or Thomson event, packet move,
   and subsequent real trace. Verify same-sign persistence, sign-change
   reset, and a new history after emission into a different line. Check event
   and estimator outcomes as scientific evidence; keep exact cursor-ID checks
   as separate orchestration assertions.
6. Test one monotonic ray with two distinct, resolvably close line frequencies.
   In the positive-gradient shell from the red case, derive frequencies from
   the unsquared `w(s)` at `s_1 = 1e12` and `s_2 = s_1 + 1e7` cm. Their
   frequency separation is much larger than floating-point spacing and the
   existing range padding, while their distances are much farther apart than
   the spatial cutoff. Assert two distinct positions and two contributions.
   This detects a cursor that skips the next physical line.
7. Exercise the zero finder with a sign-changing zero near a shell edge and a
   near-tangent pair. Use the analytic `G` and frequency extremum to decide
   whether each interval is representably resolved. A distinct forward zero
   with a frequency excursion above its propagated roundoff bound must retain
   the reversal and permit its second resonance, even when the old radial
   cutoff would omit the zero. If a zero rounds onto an endpoint or the
   spectral excursion lies within that bound, record the ray-coordinate
   separation and current solver cutoff and label the result unresolved. A
   tangent with no sign change must not reset the cursor.
8. Implement the cursor, explicit-direction range helper, endpoint rule, and
   any required zero-filter correction in the real-packet path. Add type hints
   and update touched docstrings. Run the new tests, closest nonhomologous
   modules, and Ruff. A test that already passes before implementation is
   existing-behavior coverage; do not change production code to force a diff.

The reversal positions above are separated from the zero by at least `4.6e13` cm,
so the Sobolev estimator is well conditioned there. At `nu_line` near `1e15`
Hz, one floating-point frequency step is about `0.125` Hz. The corresponding
projected-velocity uncertainty is about `c * 0.125/1e15` cm/s, which gives a
root-position scale below `100` cm at the two stated `|G|` values. Use
`atol = 100` cm for those positions, with no generic near-zero-gradient
tolerance. Compare estimator sums at about `1e-10` relative tolerance after
checking their sensitivity to that root-position uncertainty; tighten if the
implemented solver supports it. Compare event types, counts, shell IDs, and
zero additional estimator increments exactly. For the red-case boundary
distance, `rtol = 1e-12` allows roughly the floating-point error of its
geometric square root without obscuring the centimeter-scale replay event.

| Test | Claim and regime | Independent verification | Intended defect; current red status | Evidence |
| --- | --- | --- | --- | --- |
| Two-trace coherent replay | One line, positive `G` before and after scattering | Analytic gradient sign and shell-boundary distance; one line estimator contribution | Current two-trace numerical reproduction returns a second line at about 1.891 cm; test not yet added | Verification |
| Two-root reversal | One line, nonradial ray, two resolved roots and nonzero local `G` | Unsquared `w(x)` roots, analytic gradients, energy and blueward estimator equations above | A cursor retained across a real reversal; current distance solver returns roots within about 25 cm of the independent positions | Verification |
| Equal-frequency range direction | Rounded endpoint frequencies, known nonzero `G` | Explicit analytic gradient sign and ordered line IDs | Direction inferred from rounded endpoint comparison | Verification of helper contract |
| Continuous shell-edge line | Exact edge resonance, different shell opacities | Upstream ownership and one estimator contribution | Line missed or counted twice at shell edge | Verification |
| Event transitions | Coherent, Thomson, frequency-changing emission, and shell sign change | Actual event, move, next trace, and expected physical crossing sequence | Frontier lost, retained across a reversal, or retained after a frequency jump | Verification; separate ID checks are orchestration |
| Two close lines | Positive `G`, two distinct frequencies derived at known positions | Analytic `w(s)` positions and two estimator contributions | Overaggressive cursor or distance cutoff hides a physical second line | Verification |
| Near-edge and near-tangent zeros | Resolved or explicitly bounded subintervals | Analytic roots, sign change, and representable frequency excursion | Zero filter omits a resolvable reversal | Numerical verification |

For falsification, temporarily remove the cursor clamp and confirm the red
two-trace assertion fails again. Temporarily retain the cursor across the
analytic zero and confirm the second reversal crossing disappears. Revert
these mutations after each check. Search scientific papers on nonmonotonic
Sobolev transport for additional benchmark fields and resonance sequences.

## Validation and handoff

Run Python, pytest, and Ruff inside the existing `tardis` Conda environment.
In this workspace, `conda run -n tardis` reported `NoWritableEnvsDirError`,
while activating `/home/connor/micromamba/envs/tardis` through `conda.sh`
worked. Do not treat the environment as unavailable based on the failed
`conda run` invocation.

After the implementation, run the new test nodes first, then their containing
modules and Ruff on each changed Python path. Confirm the changed Numba tracer
compiles in nopython mode. Use `conda run --no-capture-output -n tardis` where
that command works, or activate the existing environment before running the
same checks in this workspace. Run the broader `pytest tardis` suite after
targeted validation unless it cannot complete in the available environment
or time budget; record the reason if omitted.

The existing `test_nonhomologous_trace_packet` and
`test_nonhomologous_trace_packet_no_line_fallthrough` cases in
`tardis/transport/montecarlo/tests/test_transport.py` and
`test_nonhomologous_packet_propagation_dispatch_numba_disabled` in
`tardis/transport/montecarlo/tests/test_single_packet_loop.py` use
`regression_data`. The dispatch test checks event routing; it is not
scientific evidence of replay prevention. Follow the repository regression
validation procedure and run these matching cases with
`--tardis-regression-data` against an external, materialized data checkout.
Do not update regression references without explicit authorization. Run the
broader TARDIS suite after targeted checks if the environment and available
time permit, and report any omitted checks.

Document the user-facing behavior with the implementation. **Do not build the
documentation**, as requested for this task.
