# Step 2: Shared Radiative Transport Core

## Summary
Expand phase 2 into a concrete refactor of radiative packet transport. The goal is to move duplicated classic/IIP/nonhomologous transport mechanics into shared jitted helpers while keeping thin mode-specific entrypoints for numba specialization and physics differences. Step 1 owns all characterization and regression testing, so this step does not include a test plan.

## Shared Modules And Interfaces
- Add shared radiative transport helpers:
  - `tardis/transport/montecarlo/path_trace.py`: shared event-candidate and line-tracing helpers.
- Add shared event types and event-candidate/result tuple helpers in `tardis/transport/montecarlo/transport_events.py`.
- Reuse existing interaction modules for shared interaction-event handling:
  - `tardis/transport/montecarlo/interaction_events.py`: shared radiative line, continuum, macro-atom, and electron-scattering dispatch helpers.
  - `tardis/transport/montecarlo/interaction_event_callers.py`: compatibility shims and mode-facing callers while imports are migrated.
- Add shared packet movement helpers directly under the existing `tardis/transport/montecarlo/packets/` package:
  - `radiative_movement.py`: shared radiative `move_r_packet` implementation.
  - `boundary_movement.py`: shared shell-boundary crossing/status handling.
  - `frame_transformations.py`: packet-level lab/CMF transforms used during movement and interactions.
- Keep existing mode entrypoints in:
  - `modes/classic/rad_packet_transport.py`
  - `modes/iip/rad_packet_transport.py`
  - `modes/nonhomologous/rad_packet_transport.py`
- Mode files should become wrappers that assemble mode-specific inputs, call shared helpers, and preserve existing numba-specialized signatures where needed.

## Planned File And Directory Structure
- Add files only inside existing directories; do not create new subdirectories for this stage.
- Add `tardis/transport/montecarlo/path_trace.py` for event candidate construction, candidate selection, and line scanner helpers.
- Add `tardis/transport/montecarlo/transport_events.py` for `TransportEvent` values and candidate/result tuple constants.
- Update existing `tardis/transport/montecarlo/interaction_events.py` and `tardis/transport/montecarlo/interaction_event_callers.py` for shared radiative line, continuum, macro-atom, and electron-scattering dispatch helpers.
- Add packet movement modules directly under `tardis/transport/montecarlo/packets/`:
  - `radiative_movement.py`: shared radiative packet position, angle, and bulk-estimator movement.
  - `boundary_movement.py`: shared shell-boundary movement and packet status updates.
  - `frame_transformations.py`: packet-level lab-to-comoving and comoving-to-lab transform helpers.
- Update existing radiative mode modules:
  - `tardis/transport/montecarlo/modes/classic/rad_packet_transport.py`: classic thin wrappers around shared trace, interaction, and movement helpers.
  - `tardis/transport/montecarlo/modes/iip/rad_packet_transport.py`: IIP thin wrappers plus continuum-specific candidate wiring.
  - `tardis/transport/montecarlo/modes/nonhomologous/rad_packet_transport.py`: nonhomologous thin wrappers plus nonhomologous line scanner wiring.
- Update `tardis/transport/frame_transformations.py` only for shared velocity-based Doppler primitives used by packet movement and interaction helpers.
- Keep `modes/*/packet_propagation.py`, `modes/*/montecarlo_transport.py`, and solver modules structurally unchanged in this stage except for imports needed to call the new shared helpers.

## Thin Numba Entrypoint Design
- Keep every public jitted mode entrypoint as a small `@njit(**njit_dict)` function with a concrete, mode-specific signature. These entrypoints should contain argument normalization, calls to shared helpers, and mode-specific dispatch only.
- Do not pass Python callables, strategy objects, dictionaries, optional callback hooks, or dynamically typed configuration into jitted helpers. Encode mode differences through separate thin entrypoints and separate shared helper calls.
- Build entrypoints in three layers:
  - mode wrapper: current import-compatible function in `modes/*/rad_packet_transport.py`;
  - shared orchestration helper: jitted helper in `path_trace.py`, `interaction_events.py`, or `packets/*.py`;
  - physics-specific primitive: homologous, nonhomologous, continuum, or interaction helper with concrete array/jitclass arguments.
- Keep return values numba-friendly: tuples of scalars, enum integers, and existing jitclass/struct values. Avoid Python dataclasses or namedtuples inside nopython code.
- Prefer duplicated two-line mode wrappers over generic abstractions that make numba infer imprecise types. The duplication target is algorithmic body code, not the existence of separate mode entry functions.
- Entrypoints should preserve current RNG draw order by sampling in the same logical layer as the old implementation. Shared helpers may sample only when every caller previously sampled at that same point.
- Each new thin entrypoint should have a direct focused test or be covered by an existing characterization test before its old implementation is removed.

## Function Moves And Modifications
- Move duplicated `move_packet_across_shell_boundary` from classic, IIP, and nonhomologous `rad_packet_transport.py` into `packets/boundary_movement.py` unchanged in behavior.
- Replace duplicated classic/IIP `move_r_packet` with one shared homologous helper:
  - New signature: `move_r_packet(r_packet, distance, geometry, estimators_bulk, enable_full_relativity)`.
  - It obtains local velocity from geometry instead of `time_explosion`.
  - Classic and IIP wrappers pass `NumbaRadial1DGeometry`.
- Modify nonhomologous `move_r_packet` to use the same shared movement body:
  - It already accepts geometry; retain nonhomologous local velocity via `geometry.get_velocity(...)`.
  - Keep nonhomologous-only velocity-gradient behavior outside the generic movement math.
- Place shared radiative movement code in `packets/radiative_movement.py`, not in `RPacket` or packet collection definition files. This keeps packet data structures separate from transport-step movement algorithms while keeping packet movement discoverable under `packets`.
- Add packet-level transform helpers in `packets/frame_transformations.py`:
  - `transform_packet_lab_to_comoving_frame(r_packet, geometry, enable_full_relativity) -> tuple[float, float, float]`.
  - `transform_packet_comoving_to_lab_frame(r_packet, geometry, comoving_nu, comoving_energy, comoving_mu, enable_full_relativity) -> None`.
  - These replace repeated inline Doppler/inverse-Doppler blocks in `thomson_scatter`, `continuum_event`, `line_scatter_event`, `line_emission`, `bound_free_emission`, `free_free_emission`, packet initialization, and opacity calculations.
- Modify `tardis/transport/frame_transformations.py`:
  - Add velocity-based Doppler helpers used by transport internals.
  - Keep old `time_explosion`-based helpers temporarily only for callers outside this refactor.
  - Route homologous transport through geometry-derived velocity, not direct `time_explosion`.

## Trace Packet Cleanup
- Split `trace_packet` into shared helpers:
  - `sample_event_optical_depth()`
  - `calculate_boundary_event_candidate(...)`
  - `calculate_continuum_or_electron_event_candidate(...)`
  - `select_nearest_trace_candidate(...)`
  - `update_packet_line_indices_after_trace(...)`
  - `scan_homologous_line_candidates(...)`
  - `scan_nonhomologous_line_candidates(...)`
- Classic and IIP should share the homologous line scanner; IIP adds continuum event selection using `chi_continuum` and `escat_prob`.
- Nonhomologous keeps its forward/backward line traversal and `calculate_distance_line_nonhomologous`, but uses the same candidate comparison and event return shape.
- Modify `calculate_distance_line` so homologous line distance is derived from geometry velocity/radius information rather than a direct `time_explosion` parameter. Keep full-relativity math behavior equivalent.
- Build trace entrypoints as:
  - `trace_packet_classic(...)`: homologous line candidates plus electron/boundary candidates;
  - `trace_packet_iip(...)`: homologous line candidates plus continuum/electron/boundary candidates;
  - `trace_packet_nonhomologous(...)`: nonhomologous line candidates plus electron/boundary candidates.
- All three trace entrypoints should return the same event tuple shape, for example `(event, distance, interaction_id, shell_change, tau_event, tau_trace)`, even if some fields are sentinel values for a mode.
- Shared candidate comparison should accept only precomputed candidate distances and event codes. It should not know how line, continuum, electron, or nonhomologous candidates were computed.

## Step-By-Step Implementation
1. Add the shared files without changing behavior:
   - Create `tardis/transport/montecarlo/transport_events.py`.
   - Create `tardis/transport/montecarlo/packets/radiative_movement.py`.
   - Create `tardis/transport/montecarlo/packets/boundary_movement.py`.
   - Create `tardis/transport/montecarlo/packets/frame_transformations.py`.
   - Create `tardis/transport/montecarlo/path_trace.py`.
   - Reuse existing `tardis/transport/montecarlo/interaction_events.py` and `tardis/transport/montecarlo/interaction_event_callers.py`; do not add a new interaction subpackage.

2. Define shared event types before moving any transport logic:
   - Move or mirror `InteractionType` semantics into `transport_events.py` as `TransportEvent`.
   - Add tuple field index constants or helper functions for the shared trace result shape.
   - Keep compatibility aliases to existing `InteractionType` values until mode wrappers and tests are updated.

3. Move shell-boundary movement first because it is behavior-identical across modes:
   - Move `move_packet_across_shell_boundary` into `packets/boundary_movement.py`.
   - Update classic, IIP, and nonhomologous `rad_packet_transport.py` to import and re-export the shared helper.
   - Keep the old function name available from each mode module while downstream imports are migrated.

4. Add velocity-based frame primitives:
   - Add velocity-based Doppler and inverse-Doppler helpers to `tardis/transport/frame_transformations.py`.
   - Keep existing `time_explosion`-based helpers for callers outside this stage.
   - Add geometry-local velocity helpers for `NumbaRadial1DGeometry` and `NumbaNonhomologousRadial1DGeometry` as needed.

5. Add packet-level frame transform helpers:
   - Implement `transform_packet_lab_to_comoving_frame(...)` in `packets/frame_transformations.py`.
   - Implement `transform_packet_comoving_to_lab_frame(...)` in `packets/frame_transformations.py`.
   - Use geometry-derived velocity and the new velocity-based Doppler primitives.
   - Do not yet replace all callers; first make the helpers available and importable.

6. Move shared radiative packet movement:
   - Move the classic/IIP `move_r_packet` body into `packets/radiative_movement.py`.
   - Change its signature to `move_r_packet(r_packet, distance, geometry, estimators_bulk, enable_full_relativity)`.
   - Use `transform_packet_lab_to_comoving_frame(...)` for CMF frequency and energy.
   - Update nonhomologous movement to use the same shared helper while preserving nonhomologous local velocity lookup.
   - Keep mode-level wrappers in classic, IIP, and nonhomologous `rad_packet_transport.py`.

7. Move interaction-event handling into `interaction_events.py` and keep compatibility callers in `interaction_event_callers.py`:
   - Move `continuum_event` as `continuum_interaction_event`.
   - Move `line_scatter_event` as `line_interaction_event`.
   - Move `macro_atom_event` as `rad_interaction_process_sampler`.
   - Move `determine_continuum_macro_activation_idx` as `select_continuum_interaction_activation_level_id`.
   - Move `determine_bf_macro_activation_idx` as `select_bound_free_interaction_activation_level_id`.
   - Move or wrap `thomson_scatter` behind the shared interaction-event interface.
   - Replace broad local names: `destination_level_idx` becomes `source_block_activation_id`; line-specific `activation_level_id` becomes `line_interaction_activation_level_id`.

8. Replace repeated interaction frame-transform code:
   - Update `continuum_interaction_event`, `line_interaction_event`, `thomson_scatter`, `line_emission`, `bound_free_emission`, and `free_free_emission` to use `packets/frame_transformations.py`.
   - Pass geometry into these functions instead of `time_explosion`.
   - Preserve RNG draw order by keeping random draws in the same logical locations as the current implementation.

9. Extract shared trace primitives into `tardis/transport/montecarlo/path_trace.py`:
   - Add `sample_event_optical_depth()`.
   - Add boundary, continuum/electron, and line event candidate helpers.
   - Add `select_nearest_trace_candidate(...)`.
   - Add `update_packet_line_indices_after_trace(...)`.
   - Add `scan_homologous_line_candidates(...)` for classic and IIP.
   - Add `scan_nonhomologous_line_candidates(...)` for nonhomologous forward/backward traversal.

10. Update line-distance and estimator dependencies used by tracing:
    - Modify homologous line-distance calculation so it receives geometry-derived radius/velocity information instead of direct `time_explosion`.
    - Keep nonhomologous `calculate_distance_line_nonhomologous` as the primitive for nonhomologous trace scanning.
    - Keep estimator updates in the same line-scan positions to preserve behavior.

11. Replace each mode `trace_packet` body with a thin wrapper:
    - Classic wrapper calls homologous line candidates plus electron and boundary candidates.
    - IIP wrapper calls homologous line candidates plus continuum/electron and boundary candidates.
    - Nonhomologous wrapper calls nonhomologous line candidates plus electron and boundary candidates.
    - All wrappers return the shared event tuple shape while preserving mode-specific signatures until callers are migrated.

12. Update packet propagation imports and calls:
    - Update classic, IIP, and nonhomologous `packet_propagation.py` to import shared event handlers, movement helpers, and trace wrappers.
    - Replace `time_explosion` arguments in transport-internal calls with geometry where possible.
    - Leave solver and `montecarlo_transport.py` structure unchanged except for any necessary argument plumbing.

13. Remove duplicate implementations once all mode wrappers delegate correctly:
    - Delete duplicated movement bodies from mode `rad_packet_transport.py`.
    - Delete duplicated trace bodies from mode `rad_packet_transport.py`.
    - Delete or reduce old `interaction_event_callers.py` and duplicate nonhomologous interaction caller modules to compatibility imports if needed.
    - Keep public/import-compatible names until all internal imports and tests use the new paths.

14. Finish the `time_explosion` boundary cleanup:
    - Remove `time_explosion` parameters from shared movement helpers, packet-level frame transforms, trace helpers, interaction handlers, and line-distance helpers.
    - Update packet initialization inside classic and IIP propagation so `set_packet_props_partial_relativity`, `set_packet_props_full_relativity`, and `RPacket.initialize_line_id` use geometry-derived velocity or packet-level frame transform helpers.
    - Update opacity calculations in packet propagation so CMF frequency comes from `transform_packet_lab_to_comoving_frame(...)`, not direct calls to `get_doppler_factor(..., time_explosion, ...)`.
    - Update virtual-packet calls only where required by changed function signatures; otherwise leave virtual-packet `time_explosion` cleanup for a later virtual-packet-specific refactor.
    - Keep `time_explosion` only where geometry, packet sources, model state, or explicitly out-of-scope virtual-packet code are constructed.
    - Leave virtual-packet and packet-source cleanup for later unless required by this stage.

## Time Explosion Boundary
- Treat `time_explosion` as an input that sets homologous radii and velocities when geometry is constructed. After `HomologousRadial1DGeometry.to_numba()` creates `NumbaRadial1DGeometry`, transport should read radius/velocity information from that geometry object.
- Add a numba-compatible homologous velocity accessor, for example `NumbaRadial1DGeometry.get_velocity(r, shell_id)`, so shared helpers can use the same interface shape as `NumbaNonhomologousRadial1DGeometry.get_velocity(...)`.
- Use geometry in these transport-internal call chains:
  - packet initialization: replace `set_packet_props_*_relativity(r_packet, time_explosion)` with geometry-based helpers;
  - line index initialization: replace `RPacket.initialize_line_id(opacity_state, time_explosion, enable_full_relativity)` with a geometry-based equivalent;
  - opacity setup: compute CMF frequency through packet-level frame transforms;
  - trace: pass geometry to line-distance and estimator-update helpers;
  - movement: pass geometry to `move_r_packet`;
  - interactions: pass geometry to line, continuum, electron-scattering, and emission helpers.
- Keep legacy `time_explosion`-based functions in `tardis/transport/frame_transformations.py` temporarily for code outside this stage and for compatibility tests. New shared radiative transport helpers should not call them.
- Do not remove `time_explosion` from solver or transport-state objects in this stage. Solvers may still receive model `time_explosion` and use it to construct geometry, configure packet sources, or preserve existing state outputs.
- Leave packet-source behavior, blackbody-relativistic source configuration, gamma-ray transport, and virtual-packet internals out of scope unless a signature change in this stage forces a narrow update.

## Assumptions
- `RPacket` remains lab-frame state; CMF values are temporary values returned by transform helpers or passed back into reverse transforms.
- Existing external imports can be updated; internal transport compatibility is more important than preserving old module paths.
- Numerical behavior and RNG draw order should remain unchanged unless step 1 characterization shows an intentional, documented difference.
