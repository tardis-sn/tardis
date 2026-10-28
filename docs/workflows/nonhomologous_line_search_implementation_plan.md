# Nonhomologous line search: traversal state implementation plan

## Science requirement and scope

A real packet must encounter a given line at most once while its projected
velocity changes monotonically along a free-flight trajectory. A line may be
encountered again after a genuine reversal of the projected velocity gradient.
Count a resonance as traversed once its Sobolev optical depth and line
estimators have been processed, whether or not it causes a scattering event.

This plan has been revised after auditing implementation commit `d5d61544e8`
against its parent, `b0bf479057`, which committed the previous plan following
[the static audit](nonhomologous_line_search_plan_audit.md). The audit changes
were followed by a scoped implementation simplification: the cursor remains,
while the new extremum exception band, edge-ownership policy, and zero-finder
refinements have been removed. The distance helpers are restored to their
pre-fix implementation, including their existing direction calculation.
The findings below cite the original commit;
the separately scoped concerns explain which requirements of the previous
plan have been withdrawn. The core fix covers the line search in
`tardis/transport/montecarlo/modes/nonhomologous/rad_packet_transport.py`.
The packet loop in `modes/nonhomologous/packet_propagation.py` must carry search
state between calls to `trace_packet`. The separate nonhomologous virtual-packet
tracer uses the sign of radial `dv/dr` and does not split at projected-gradient
zeros. Assess that path separately. Until then, this fix establishes
real-packet line-search behavior; it does not establish corrected spectra.

Assumptions: shell velocity is piecewise linear in radius; the transport line
list is sorted by descending rest frequency; nonhomologous real-packet
transport uses its existing first-order Doppler treatment. The Sobolev formula
is singular at exactly zero projected gradient; robust transport through that
isolated coincidence is outside this fix. No issue identifier was supplied.

The user's scope clarification accepts losing an isolated packet, with its
luminosity contribution, when its angle, location, or frequency must coincide
with a degeneracy to machine precision. Prefer terminating that packet and
continuing over adding general machinery or aborting the simulation for such
a case. No new packet-loss framework is introduced by this simplification;
pre-existing singular-resonance handling is retained. Ordinary finite-gradient
resonances still require correct processing. No packet-loss rate or luminosity
error has been measured.

## Audit findings: why the implementation grew

The diagnosis remains simple: each trace rebuilds a padded frequency range
without remembering which line was processed. An identity-based cursor fixes
that defect; a larger distance epsilon cannot reliably do so. The quartic
solver, binary searches, range padding, and analytic interval splitting
already existed before `d5d61544e8` and are not new complexity from this attempt.

Scope growth occurred at two stages. The preceding plan already required
exact edge ownership, near-tangent zero recovery, propagated roundoff analysis,
and an extensive event matrix. The implementation mostly followed those
expanded requirements. It then went further by introducing a new runtime
exception region around frequency extrema. Of its 1,723 inserted lines,
1,493 are tests; test volume and production complexity should be assessed
separately.

The principal findings, with locations in `d5d61544e8`, are:

1. **Remove the new near-extremum abort policy from the replay fix.**
   `modes/nonhomologous/rad_packet_transport.py:136–160,209–217,267–271`
   computes a frequency uncertainty band, assigns every candidate in that
   band the extremum's distance without solving for a resonance, and raises.
   The new `test_nonhomologous_line_below_projected_frequency_minimum_is_bounded`
   independently proves with Decimal arithmetic that its near line has no
   physical resonance, then requires an exception anyway
   (`tests/test_nonhomologous.py:657–745`). This validates a new failure policy,
   not a necessary consequence of preventing replay. The previous plan asked
   to label unresolved evidence and preserve the exact-zero singularity;
   neither requirement authorizes inventing a resonance for a nearby line.
   The clarified scope also rejects making a machine-precision coincidence
   a simulation-wide failure condition.
2. **The uncertainty threshold is not a demonstrated error bound.**
   `gamma_32` assumes an operation count without deriving an error budget for
   the calculation, geometry, or extremum location. Its leading term at
   `nu_lab = 1e15` Hz is about `3.55` Hz, or 28 frequency ULPs, before the
   velocity terms. Tests at zero/one ULP and at an offset of `-8` Hz do not
   establish this threshold. Passing assertions that demand the exception
   cannot supply its missing scientific justification.
3. **Geometric resolution was promoted into a transport requirement.**
   The near-start and near-tangent tests
   (`tests/test_nonhomologous.py:1186–1363`) explicitly limit their claims to
   zero locations. The tangent example's frequency excursion is below
   `1e-9` Hz, versus `0.125` Hz spacing. It cannot demonstrate a resolvable
   second line encounter. These tests do not justify coupling the Newton
   refinement of `cbrt`, the zero finder's changed fastmath policy, or a runtime
   rejection region to the replay fix. A resolved transport failure would be
   needed to bring those numerical changes into scope.
4. **An obsolete branch remains after the boundary tie change.**
   In `modes/nonhomologous/rad_packet_transport.py:254–266`, the check
   `distance_boundary < distance_trace` is unreachable: the earlier check
   rejects `distance_trace > interval_end`, and every interval ends at or
   before the boundary. The old `<=` branch handled equality; changing that
   rule left dead code behind.

The two scalar state values are not the source of needless abstraction. Their
event lifecycle is necessary once the fix remembers prior crossings: genuine
gradient reversals and frequency-changing emission must allow a fresh search.
The removed midpoint-to-quarter-point fallback had a mathematical case: a
tangent with no sign change can place exactly zero gradient at the midpoint.
That finely tuned sample is outside the robustness requirement for this fix;
it did not demonstrate an additional failure for ordinary packets.

Continuing the scope review during simplification removed two more layers.
The explicit-direction range-helper API and interval-gradient sampling were
only needed to extend direction handling to equal or nearly equal rounded
endpoint frequencies. Existing endpoint-derived direction suffices for the
resolved monotonic intervals in scope. After a line event, unconditional
reseeding from the emitted line and outgoing gradient also removes the
absorbed-line comparison: for coherent same-sign emission it assigns exactly
the cursor value already returned by the tracer.

All abbreviated production and test paths above are relative to
`tardis/transport/montecarlo/`.

### Physical relevance and priorities

Exact coincidences are constructible in a test, but that does not make them
representative of sampled packets. Conversely, an event imposed by the
algorithm must not be dismissed on probabilistic grounds: line emission
deliberately sets a packet's comoving frequency to a line frequency. The
review therefore uses the following priorities, rather than trying to cover
every combination of floating-point inputs.

| Situation | Relevance to actual packets | Scope decision |
| --- | --- | --- |
| A just-processed or emitted line appears a few centimeters ahead on the next trace | Frequency equality is imposed by emission, and the solver can perturb the distance; this is the reported bug | Fix and retain a real packet-loop regression |
| A line reappears after an ordinary projected-gradient reversal | A finite interval of rays and line frequencies can have two well-separated, finite-gradient roots | Keep one well-conditioned analytic example |
| The next distinct line is close to the previous line | Dense line lists naturally contain nearby transitions | Keep one resolvable close-line example; no arbitrary-separation guarantee |
| Coherent or Thomson scattering changes the ray; a different line is emitted | Normal event types change the validity of the cursor | Keep their simple state transitions and focused evidence |
| A packet crosses an ordinary shell edge | Every escaping packet crosses mesh boundaries; retaining processed identity is relevant | Preserve history for continuous frequency and matching direction; reset simply for a velocity jump |
| A resonance lies exactly on a shell edge | The boundary crossing is common, but equality of a particular resonance and edge is a separate finely tuned condition | Do not require a new ownership convention or residual-clamping infrastructure for this fix |
| A resonance lies exactly at `G = 0`, especially also exactly on an edge | Coincidence of the line frequency with the trajectory's extremum is finely tuned; adding edge equality is still more restrictive | No robust-resolution requirement; isolated packet loss is acceptable |
| Geometric zero pairs or endpoint intervals have spectral excursion below representable frequency resolution | Their existence in ray coordinates does not establish a resolvable physical line encounter | Document the limit; do not refine geometry solely to make these synthetic cases pass |

Finite-gradient resonances near a turning point are not all irrelevant. Their
conditioning can matter, and genuine reversals must work. What is deferred is
exhaustive treatment of machine-precision degeneracies, not the nonmonotonic
transport problem as a whole. This audit does not estimate occurrence rates
from production simulations and does not assume that a broad frequency band
is equivalent to an exact coincidence.

## Reproduced defect

For one shell with radii `[1e14, 2e14]` cm and edge velocities
`[1e8, 1.8e8]` cm/s, take a packet at `r = 1.25e14` cm, `mu = 0.1`, and a
line at `1e15` Hz. Set the lab-frame packet frequency to
`nu_line / (1 - v(r) * mu / c)`. The field is genuinely nonhomologous:
`v(r) = (8e-7 s^-1) r + 2e7 cm/s`. Its projected gradient stays positive
along this ray, so there is no physical opportunity to revisit the line.

As a direct diagnostic, set `next_line_id = 1` and `prev_line_id = 0` to
indicate that this line was crossed. The pre-fix `trace_packet` still selects it
about **1.455 cm** ahead and adds a line-estimator contribution. This packet-ID
setup is not the primary red test: those IDs will not be the new traversal
state.

The primary red case uses two actual traces. In the same shell, launch from
`r = 1.25e14` cm and `mu = 0.1`. Construct the lab frequency so the first
resonance of the `1e15` Hz line lies `1e12` cm along the ray. Give it high line
optical depth and negligible electron opacity. After the first `LINE` event,
move the packet to that resonance, set the coherent outgoing direction to
`mu = 1`, and emit at the same line. On the pre-fix implementation the second
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

The zero finder has numerical resolution limits near tangency and interval
endpoints. These are valid limitations to document, but resolving every
representable geometric reversal is not an acceptance condition for the replay
fix. The added geometry tests do not establish a lost, spectrally resolved
resonance. Keep the ordinary two-root physics test well away from `G = 0`;
assess zero-filter and root-conditioning improvements separately, as described
under the vetted concerns below.

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

1. Reuse `get_line_id_range_nonhomologous` and its returned `line_id_step`.
   For a resolved monotonic interval, the difference between endpoint
   frequencies gives the traversal direction. No new helper API or interval
   sampling is needed.
2. If the active direction matches `line_search_direction`, intersect the
   directional ID range with the cursor. For redshift, start at least at
   `next_untraversed_line_id`; for blueshift, start at most there. If direction
   is `0` or has reversed, start from the frequency-derived range.
3. Solve distances only for eligible lines. Advance the cursor to
   `line_id + direction` immediately after processing a real resonance, before
   any return for a line event. Do not mark a candidate if electron scattering
   or a boundary occurs first.
4. When an analytic zero reverses `G`, clear the old frontier and derive a new
   range for the next interval. Do not synthesize an event there for every
   candidate near the extremum. Exact singular coincidences have the limited
   scope and accepted packet-loss fallback stated above.

Equal or nearly equal rounded endpoint frequencies retain the pre-existing
helper's numerical limitations. Extending that resolution is outside this fix.

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

After a line event, seed the outgoing cursor from the emitted line ID set by
`line_emission` and the outgoing gradient sign. For coherent same-sign emission
this equals the already advanced cursor, so no absorbed-line comparison is
necessary. Update after `line_scatter_event` has chosen the new direction and
frequency, before calling `trace_packet` again. After Thomson scattering,
the next trace can compare its outgoing interval direction with the stored
direction. A separate eager reset is unnecessary if that comparison already
handles a reversal.

### Endpoint scope

The previous plan required a new upstream-shell ownership convention for exact
edge resonances. The audit identifies this as a separate, pre-existing missing
resonance concern, rather than a prerequisite for preventing replay. Defer that
convention and its solver changes from the narrow fix. Cursor persistence
across an ordinary continuous shell crossing remains in scope. This plan does
not claim complete endpoint coverage or improved root resolution.

Keep the existing distance tolerance for numerical root filtering. The cursor
provides the identity-based guard against replay; increasing a spatial
epsilon could suppress a physically distinct close line. The current solver
already rejects roots inside its `1e-14 * shell_width` spatial cutoff, so this
change does not promise resolution of arbitrarily close physical lines.

## Retained implementation and validation plan

The simplified implementation retains the following physical claims and
tests. Execution and mutation checks below describe subsequent validation;
they are not instructions to run tests under the user's current restriction.

1. Retain `test_nonhomologous_packet_loop_retains_coherent_line_cursor` as the
   primary regression, with the real tracer and real packet move.
   Use the two-trace scenario above, a controlled `mu = 1` coherent outgoing
   direction, one high-opacity line, and negligible electron opacity. Fix the
   angle draw at its sampling boundary or use an equivalent deterministic
   fixture; do not mock `trace_packet`. Confirm that the pre-fix code fails
   because it records a second `LINE` event or estimator contribution after
   the first move and emission. The implemented path must instead reach the
   analytically calculated outer boundary with exactly one line contribution.
   Bound the test's event recording so the pre-change replay cannot spin
   indefinitely; record real tracer results rather than replacing them.
   The existing `test_nonhomologous_packet_propagation_dispatch_numba_disabled`
   mocks transport dispatch and cannot establish this claim. A direct tracer
   test that injects explicit cursor state may supplement the packet-loop test
   after the cursor contract exists; an API mismatch is not a valid red phase.
2. Retain `test_nonhomologous_line_reappears_after_gradient_reversal` with the
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
3. Reuse the existing range helper without extending its API or adding
   machine-precision endpoint-direction tests.
4. Test same-sign preservation using exactly shared velocity edges, and test
   sign-change reset separately. For an allowed velocity jump, assert only
   that the cursor clears; do not claim a continuous-resonance result.
5. Retain the physical line-emission test and its separate packet-loop state
   transfer check. Same-sign Thomson cursor preservation is checked statically
   below; it does not need a second large integration scenario. Keep exact
   cursor-ID checks separate from physical estimator evidence.
6. Test one monotonic ray with two distinct, resolvably close line frequencies.
   In the positive-gradient shell from the red case, derive frequencies from
   the unsquared `w(s)` at `s_1 = 1e12` and `s_2 = s_1 + 1e7` cm. Their
   frequency separation is much larger than floating-point spacing and the
   existing range padding, while their distances are much farther apart than
   the spatial cutoff. Assert two distinct positions and two contributions.
   This detects a cursor that skips the next physical line.
7. Keep endpoint ownership and sub-resolution zero-finder studies separate
   from these acceptance tests. A tangent with no sign change must not reset
   the cursor, but recovering spectrally unresolved zero pairs is not a
   requirement of this fix.
8. Keep the cursor in the real-packet path with the event lifecycle above.
   Add type hints
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

| Existing test in `d5d61544e8` | Claim and regime | Independent verification | Intended defect and scope decision | Evidence |
| --- | --- | --- | --- | --- |
| `test_nonhomologous_packet_loop_retains_coherent_line_cursor` | One line, positive `G` before and after coherent scattering | Real tracer, move, emission, analytic boundary and one contribution | Lost history replays a line; keep as the primary regression | Verification |
| `test_nonhomologous_coherent_line_does_not_replay` | Same ray in two direct traces | Analytic first position and boundary; unchanged deposition after second trace | Removed duplicate scenario; the actual-loop regression covers replay | Verification |
| `test_nonhomologous_line_reappears_after_gradient_reversal` | One line, two well-separated roots with finite local `G` | Unsquared `w(x)` roots, analytic gradients and estimator equations above | Cursor incorrectly retained across a reversal; keep | Verification |
| `test_nonhomologous_shell_gradient_reversal_revisits_line` | Two velocity-continuous shells with opposite radial gradients | Both profiles reach the same projected velocity at independently prescribed radii | Cursor incorrectly retained when entering a new shell with reversed drift; keep | Verification |
| `test_nonhomologous_search_resolves_two_close_lines` | Positive `G`, two distinct resolvable frequencies | Analytic positions and two estimator contributions | Off-by-one cursor hides the next physical line; keep | Verification |
| `test_nonhomologous_thomson_scatter_keeps_same_direction_line_frontier` | Actual elastic event between traces, same gradient sign | Prescribed positions, Doppler energy relation and contributions | Removed after PR review; existing Thomson/dispatch tests plus the static state-transfer check below are sufficient for this fix | Previously verification and orchestration |
| `test_nonhomologous_frequency_changing_emission_starts_new_line_history` and `test_nonhomologous_packet_loop_reseeds_after_different_line_emission` | Different emitted line, new frequency history | Analytic positions with manually reseeded cursor; separate mocked packet-loop transfer | Old history hides a physical crossing; useful layered coverage, not end-to-end macro-atom verification | Verification plus orchestration |
| `test_nonhomologous_range_uses_projected_gradient_direction` | Equal or one-ULP-separated endpoints, supplied direction | Known integer order of descending line IDs | Removed with the precision-case helper API extension | Helper-contract verification |
| Shell-edge ownership tests in `test_transport.py` | Exact resonance/edge coincidence with different shell opacities | Chosen upstream convention and estimator channels | Removed from replay fix with the new ownership convention | Convention verification |
| Near-start and near-tangent zero-finder tests | Geometric reversals with unresolved spectral excursions | Analytic cubic, gradient signs; `1 cm` / `2e5 cm` position tolerances | Removed; the zero finder is restored to its pre-fix implementation | Numerical characterization |
| Exact/one-ULP extremum and below-minimum exception tests | Fine-tuned frequencies at a trajectory's extremum | Decimal establishes the below-minimum line is unreachable; exception assertions encode new policy | Removed with the unsupported exception-band policy | Characterization plus policy assertions |

The tolerances for the retained physical tests are justified by the
well-conditioned positions and estimators above; exact line-ID assertions
check orchestration only. Red-phase and mutation/falsification evidence was
not established in this audit. The recorded pre-fix centimeter distances are
the original plan's diagnostic observations, not new audit measurements.

Coverage also remains narrower than the previous plan requested: the commit
does not contain physical sign-reversing Thomson/coherent event cases or a
velocity-jump packet-loop reset case. Those are relevant cursor-lifecycle
checks; prioritize a focused example if later evidence warrants it, rather
than an exhaustive event matrix. The removed rounded-frequency test only
established integer ordering, not a resolved physical crossing.

For falsification, temporarily remove the cursor clamp and confirm the red
two-trace assertion fails again. Temporarily retain the cursor across the
analytic zero and confirm the second reversal crossing disappears. Revert
these mutations after each check when implementation validation is authorized.
No mutations are part of this audit. Search scientific papers on nonmonotonic
Sobolev transport for additional benchmark fields and resonance sequences.

## Vetted additional concerns, outside the narrow fix

**Exact shell-edge resonance ownership.** The parent tracer gives a boundary
event priority over an equal-distance line event, while the distance solver
excludes the next interval's starting point. Thus an exact edge resonance can
be processed by neither shell. This is a real mathematical limitation, but
the equality is finely tuned and the user does not require robust support
here. The upstream convention in `d5d61544e8` is a new observable choice when
shell opacities or gradients differ, not something implied by replay
prevention. It also uses exact equality of rounded frequencies to bypass the
solver (`rad_packet_transport.py:218–224`), rather than the unsquared-equation
residual criterion demanded by the preceding plan. Record the limitation;
defer both the convention and extra root-admission machinery.

**Zero filtering and root conditioning.** The previous zero finder uses a
radial tolerance near tangency and a distance tolerance near interval ends.
Near tangency, `x_zero = sqrt(r_zero^2 - b^2)` amplifies a radial perturbation,
so that tolerance is not a bound on ray-coordinate separation. This is a valid
geometric observation. However, the added cases do not establish a
spectrally resolved missing line: the near-start test explicitly disclaims
frequency resolution, and the tangent test's excursion is far below one ULP.
Similarly, a resonance near a frequency extremum is ill conditioned because
`d nu_comov/ds` approaches zero. That motivates a documented numerical limit,
not comprehensive recovery of constructed degeneracies. Neither an unresolved
geometric interval nor a candidate within range padding is proof of a
physical resonance. Retain the user-accepted isolated packet-loss option;
do not make recovery of these cases or a new simulation-wide exception band
part of the replay fix.

**Shell continuity and emission history.** These are relevant qualifications
of the cursor, rather than separate solver projects. Adjacent geometry edges
may be only approximately equal (`Radial1DGeometry` uses `np.allclose`), so a
same-sign gradient alone does not establish continuous comoving frequency.
A simple cursor reset for a stored velocity jump is appropriate. Emission
into a different line likewise breaks the old contiguous frequency history.
Keep these narrow state actions; no shock-transport validation or general
invalid-input defense is needed.

## Validation and handoff

The audit and initial simplification each used three subagents for
production, tests, and independent scope review. At that stage the production
change relative to `b0bf479057` is 77 inserted and 6 deleted lines in the tracer
and packet loop. `calculate_distances.py` is identical to that parent, verified
with `git diff --exit-code`. Seven new focused behavior tests remained; 722
test lines and the superseded numerical policies were removed from the attempted
implementation. Existing untracked work, dependencies, CI, scientific
assumptions, and regression references were preserved.

### PR #3802 comment follow-up

The [review comments](https://github.com/tardis-sn/tardis/pull/3802)
requested clearer direction and cursor explanations, fixture-based setup,
plainer test wording, and reconsideration of the heavyweight Thomson test.
The follow-up adds comments explaining redshift/blueshift and the cursor
clamp, reuses scenario fixtures, and preserves the retained physical inputs,
analytic expectations, and tolerances. Six focused new tests remain after
removing the additional Thomson integration scenario.

The Thomson decision is based on this one-time static consistency check:

- Both tracer returns for electron scattering carry the current cursor and
  direction without assigning a new crossing.
- The packet loop's electron-scattering branch does not modify either scalar.
- `thomson_scatter` preserves the comoving frequency and cannot access the
  loop-owned cursor.
- The next trace keeps the frontier for matching drift direction and resets
  it when the direction reverses.

Existing Thomson interaction and packet-loop dispatch tests remain. This
checks the cursor wiring relevant to this fix; it does not claim a new
validation of all Thomson-scattering physics. No new transport behavior is
introduced by the comment follow-up, and no tests or documentation builds
were run for it.

Follow-up static validation: Ruff on the two production modules, three test
modules, and their shared `conftest.py` reported the same 88 legacy diagnostics
as `1a72ddf3b6`, with no new diagnostics. Comparing parsed syntax trees confirms
both production modules changed only in comments. `git diff --check` passed.

The user instructed that no documentation build or test suite be run for this
audit. That restriction was retained during implementation: no tests ran after
the instruction, and no documentation build was
started. Before the instruction arrived, the test reviewer had already
completed the following focused commands in the existing environment:

```shell
source /home/connor/miniforge3/etc/profile.d/conda.sh
conda activate /home/connor/micromamba/envs/tardis
NUMBA_CACHE_DIR=/tmp/tardis-review-test-audit pytest tardis/transport/montecarlo/tests/test_nonhomologous.py -k 'line_reappears or gradient_zero_is or below_projected or search_resolves or thomson_scatter_keeps or frequency_changing or range_uses or zero_finder' -q
NUMBA_CACHE_DIR=/tmp/tardis-review-test-audit pytest tardis/transport/montecarlo/tests/test_transport.py::test_nonhomologous_coherent_line_does_not_replay tardis/transport/montecarlo/tests/test_single_packet_loop.py::test_nonhomologous_packet_loop_retains_coherent_line_cursor -q
```

Results were `13 passed, 57 deselected` and `2 passed`, respectively. These are
historical results for the attempted implementation, not validation of the
simplified source or the scientific necessity of the removed abort policy.
No pre-fix red run, mutation/falsification run, complete affected
module, broader suite, or regression-data comparison was performed in this
audit. Those checks and documentation builds are omitted under the user's
instruction. Runtime and Numba compilation validation of the simplified source
remain incomplete.

Ruff was run in the activated `tardis` environment on all six affected Python
paths (the geometry module, tracer, packet loop, and three test modules).
Whole-file lint reported 70 existing annotation/docstring errors. A Ruff JSON
comparison of each file against `b0bf479057`, using the same configuration and
matching rule, message, and source line, found no new diagnostics; the parent
has 72. No unrelated annotation cleanup was added. `git diff --check` passed.
Environment checks confirmed `tardisbase`, Ruff,
and Sphinx are installed; `conda run -n tardis` failed with
`NoWritableEnvsDirError`, but explicit environment activation worked.

For subsequent implementation validation, when authorized, use the targeted
retained tests first, then affected modules, Ruff, and a Numba nopython check.
Apply the repository's broader-suite and regression-validation guidance;
neither is claimed complete by this audit.

The existing `test_nonhomologous_trace_packet` and
`test_nonhomologous_trace_packet_no_line_fallthrough` cases in
`tardis/transport/montecarlo/tests/test_transport.py` and
`test_nonhomologous_packet_propagation_dispatch_numba_disabled` in
`tardis/transport/montecarlo/tests/test_single_packet_loop.py` use
`regression_data`. The dispatch test checks event routing; it is not
scientific evidence of replay prevention. Subsequent validation should use
the repository regression-validation procedure and these matching cases with
`--tardis-regression-data` against an external, materialized data checkout.
Do not update regression references without explicit authorization.

The propagation walkthrough now documents the narrowed behavior and numerical
limits. **Do not build documentation or run further tests under the current
user restriction.**
