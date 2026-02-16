# High Energy Transport Restructuring Plan

## Context
Following PR #3427 which introduced a modular architecture for Monte Carlo transport modes (classic/iip), the High Energy (HE) transport code needs similar restructuring to align with the new architecture.

## Current State

### High Energy Transport Code Locations:
```
tardis/energy_input/
├── main_gamma_ray_loop.py           # Main loop orchestration
├── gamma_ray_transport.py           # Transport calculations
├── gamma_ray_channel.py             # Decay channel handling
├── gamma_ray_estimators.py          # HE estimators
├── transport/
│   ├── gamma_packet_loop.py         # Packet propagation loop
│   ├── gamma_ray_interactions.py    # Compton, pair creation, etc.
│   ├── gamma_ray_grid.py            # Grid/geometry for HE
│   └── GXPacket.py                  # HE packet class
└── ...

tardis/transport/montecarlo/packet_source/
└── high_energy.py                   # GammaRayPacketSource

tardis/workflows/high_energy/
└── tardis_he_workflow.py            # HE workflow
```

### New MC Transport Structure (from PR #3427):
```
tardis/transport/montecarlo/modes/
├── __init__.py
├── classic/
│   ├── __init__.py
│   ├── montecarlo_transport.py      # Main transport loop
│   ├── packet_propagation.py        # Single packet logic
│   ├── rad_packet_transport.py      # R-packet transport
│   └── solver.py                    # Mode-specific solver
└── iip/
    ├── __init__.py
    ├── montecarlo_transport.py
    ├── packet_propagation.py
    ├── rad_packet_transport.py
    └── solver.py
```

## Proposed Restructure

### Phase 1: Create HE Mode Directory Structure
Create a new `he` mode under `tardis/transport/montecarlo/modes/`:

```
tardis/transport/montecarlo/modes/he/
├── __init__.py                      # Module exports
├── solver.py                        # HE transport solver (entry point)
├── gamma_transport.py               # Main HE transport loop (from main_gamma_ray_loop.py)
├── packet_propagation.py            # Single gamma packet propagation (from gamma_packet_loop.py)
├── gamma_interactions.py            # Physics interactions (from gamma_ray_interactions.py)
├── gamma_grid.py                    # Geometry/grid handling (from gamma_ray_grid.py)
├── estimators.py                    # HE-specific estimators (from gamma_ray_estimators.py)
└── tests/
    ├── __init__.py
    └── ...                          # Migrate existing HE transport tests
```

### Phase 2: Consolidate Packet Source
Move and integrate packet source:

```
tardis/transport/montecarlo/modes/he/
└── packet_source.py                 # Move from montecarlo/packet_source/high_energy.py
```

**Changes:**
- Move `GammaRayPacketSource` from `tardis/transport/montecarlo/packet_source/high_energy.py`
- Integrate with mode-specific structure
- Update imports in dependent code

### Phase 3: Migrate Core HE Logic
Restructure files from `tardis/energy_input/` into the new mode:

#### 3.1 Create `modes/he/gamma_transport.py`
**Source:** `energy_input/main_gamma_ray_loop.py`
- Move `run_gamma_ray_loop()` function
- Move helper functions like `get_effective_time_array()`, `calculate_electron_number_density()`
- Update imports to reference new locations

#### 3.2 Create `modes/he/packet_propagation.py`
**Source:** `energy_input/transport/gamma_packet_loop.py`
- Move `gamma_packet_loop()` function
- Keep numba decorators and optimizations

#### 3.3 Create `modes/he/gamma_interactions.py`
**Source:** `energy_input/transport/gamma_ray_interactions.py`
- Move interaction functions:
  - `compton_scatter()`
  - `pair_creation_packet()`
  - `scatter_type()`
  - `get_compton_fraction_artis()`

#### 3.4 Create `modes/he/gamma_grid.py`
**Source:** `energy_input/transport/gamma_ray_grid.py`
- Move grid/geometry functions:
  - `distance_trace()`
  - `move_packet()`

#### 3.5 Create `modes/he/estimators.py`
**Source:** `energy_input/gamma_ray_estimators.py`
- Move HE-specific estimator classes/functions
- Align with new estimator architecture from PR #3427

#### 3.6 Migrate GXPacket
**Source:** `energy_input/transport/GXPacket.py`
- Move to `modes/he/packets.py` or keep in shared location
- Consider: should this stay in `energy_input/transport/` or move to modes?
- Update all imports

### Phase 4: Create HE Mode Solver
Create `modes/he/solver.py` following pattern from classic/iip modes:

```python
# modes/he/solver.py
class HighEnergyTransportSolver:
    """
    High Energy (gamma-ray) transport solver.
    
    Handles radioactive decay gamma-ray and positron transport
    through the supernova ejecta.
    """
    
    def __init__(self, ...):
        # Initialize HE-specific state
        pass
    
    def run(self, ...):
        """
        Run high energy transport simulation.
        
        Returns
        -------
        HETransportState with gamma-ray estimators
        """
        # Call gamma_transport main loop
        pass
```

### Phase 5: Update Workflows
Update `tardis/workflows/high_energy/tardis_he_workflow.py`:

**Changes:**
- Import from new location: `from tardis.transport.montecarlo.modes.he import ...`
- Update calls to use new solver interface
- Ensure backward compatibility where needed

### Phase 6: Update Imports Across Codebase
Files that will need import updates:
- `tardis/workflows/high_energy/tardis_he_workflow.py`
- `tardis/workflows/high_energy/tests/test_tardis_he_workflow.py`
- Any simulation runners that use HE transport
- Documentation notebooks

### Phase 7: Deprecation and Cleanup
After migration is complete and tested:

**Option A: Remove old locations**
```
# Mark for removal after transition period:
tardis/energy_input/main_gamma_ray_loop.py
tardis/energy_input/gamma_ray_transport.py
tardis/energy_input/transport/gamma_packet_loop.py
tardis/energy_input/transport/gamma_ray_interactions.py
tardis/energy_input/transport/gamma_ray_grid.py
```

**Option B: Keep compatibility layer**
- Keep old files as thin wrappers that import from new locations
- Add deprecation warnings
- Remove in future major version

**Recommended:** Keep `tardis/energy_input/` for decay physics and data:
- `gamma_ray_channel.py` - decay channel physics
- `decay_radiation.py` - radioactive decay data
- `energy_source.py` - energy source abstractions
- `nuclear_energy_source.py` - nuclear energy calculations
- `samplers.py` - statistical samplers
- `util.py` - shared utilities

Move transport-specific code to `modes/he/`.

## File Structure After Restructure

```
tardis/
├── energy_input/                    # Keep: decay physics & data
│   ├── gamma_ray_channel.py        # Keep
│   ├── decay_radiation.py          # Keep
│   ├── energy_source.py            # Keep
│   ├── nuclear_energy_source.py    # Keep
│   ├── samplers.py                 # Keep
│   ├── util.py                     # Keep (or move HE-specific utils to modes/he/)
│   └── transport/                  # Deprecate/remove
│       └── GXPacket.py             # Decision: move or keep?
│
├── transport/montecarlo/
│   ├── modes/
│   │   ├── classic/                # Existing (PR #3427)
│   │   ├── iip/                    # Existing (PR #3427)
│   │   └── he/                     # NEW: High Energy mode
│   │       ├── __init__.py
│   │       ├── solver.py           # Entry point
│   │       ├── gamma_transport.py  # Main loop
│   │       ├── packet_propagation.py
│   │       ├── gamma_interactions.py
│   │       ├── gamma_grid.py
│   │       ├── estimators.py
│   │       ├── packet_source.py    # From packet_source/high_energy.py
│   │       ├── packets.py          # GXPacket and related
│   │       └── tests/
│   │           └── ...
│   │
│   └── packet_source/
│       └── high_energy.py          # Remove or deprecate
│
└── workflows/
    └── high_energy/
        └── tardis_he_workflow.py   # Update imports
```

## Implementation Checklist

### Phase 1: Setup
- [ ] Create `tardis/transport/montecarlo/modes/he/` directory
- [ ] Create `__init__.py` with module documentation
- [ ] Create `tests/` subdirectory

### Phase 2: Core Transport Files
- [ ] Create `modes/he/solver.py` (HE solver class)
- [ ] Create `modes/he/gamma_transport.py` (from `main_gamma_ray_loop.py`)
- [ ] Create `modes/he/packet_propagation.py` (from `gamma_packet_loop.py`)
- [ ] Create `modes/he/gamma_interactions.py` (from `gamma_ray_interactions.py`)
- [ ] Create `modes/he/gamma_grid.py` (from `gamma_ray_grid.py`)
- [ ] Create `modes/he/estimators.py` (from `gamma_ray_estimators.py`)

### Phase 3: Packet Handling
- [ ] Create `modes/he/packet_source.py` (from `packet_source/high_energy.py`)
- [ ] Create `modes/he/packets.py` or decide on GXPacket location
- [ ] Update all GXPacket imports

### Phase 4: Test Migration
- [ ] Migrate tests from `energy_input/transport/tests/` to `modes/he/tests/`
- [ ] Migrate tests from `energy_input/tests/` (HE-specific ones)
- [ ] Update test imports and fixtures
- [ ] Ensure all tests pass

### Phase 5: Update Workflows
- [ ] Update `workflows/high_energy/tardis_he_workflow.py` imports
- [ ] Update workflow tests
- [ ] Test end-to-end HE workflow

### Phase 6: Update Documentation
- [ ] Update import paths in docstrings
- [ ] Update any HE transport documentation
- [ ] Add conceptual changelog entry (like `2026_monte_carlo_mode_architecture.md`)
- [ ] Update notebooks that use HE transport

### Phase 7: Deprecation
- [ ] Add deprecation warnings to old import locations
- [ ] Update `energy_input/__init__.py` to warn about moved functions
- [ ] Document migration path for users

### Phase 8: Final Cleanup
- [ ] Run full test suite
- [ ] Update CI/CD if needed
- [ ] Code review and approval
- [ ] Merge after PR #3427 is merged

## Benefits of Restructure

1. **Consistency**: HE transport follows same pattern as classic/iip modes
2. **Modularity**: Clear separation between decay physics (`energy_input/`) and transport (`modes/he/`)
3. **Maintainability**: All HE transport code in one location
4. **Extensibility**: Easy to add new HE transport modes or variations
5. **Testing**: Isolated test suites for each mode
6. **Documentation**: Clearer organization for users and developers

## Risks and Considerations

1. **Backward Compatibility**: Need deprecation period for old imports
2. **GXPacket Location**: Decide if GXPacket stays in `energy_input/transport/` or moves to `modes/he/`
3. **Shared Utilities**: Some utility functions might be used by both decay physics and transport
4. **Dependencies**: Ensure circular imports are avoided
5. **Testing Coverage**: Need comprehensive tests to verify migration correctness
6. **Timing**: Should be done after PR #3427 merges to avoid conflicts

## Open Questions

1. **GXPacket**: Should `GXPacket` stay in `energy_input/transport/` or move to `modes/he/packets.py`?
   - **Recommendation**: Move to `modes/he/packets.py` for consistency with other packet types

2. **Utilities**: Should HE-specific utilities in `energy_input/util.py` move to `modes/he/util.py`?
   - **Recommendation**: Move HE transport-specific utilities, keep decay physics utilities

3. **Estimators**: Should HE estimators follow the new estimator architecture from PR #3427?
   - **Recommendation**: Yes, align with `EstimatorsBulk`, `EstimatorsLine`, `EstimatorsContinuum` pattern

4. **Solver Interface**: Should `HighEnergyTransportSolver` have the same interface as `MonteCarloTransportSolver`?
   - **Recommendation**: Similar interface but adapted for HE-specific needs (different inputs/outputs)

5. **Decay Physics**: What stays in `energy_input/`?
   - **Recommendation**: Keep decay data, isotope handling, energy source abstractions
   - Move: Transport loops, packet propagation, interactions, geometry

## Timeline Estimate

- **Phase 1-2**: 2-3 days (setup and initial file creation)
- **Phase 3**: 5-7 days (migrate core logic, test)
- **Phase 4**: 2-3 days (create solver, integrate)
- **Phase 5**: 2-3 days (update workflows and tests)
- **Phase 6-7**: 2-3 days (documentation and deprecation)
- **Phase 8**: 2-3 days (final testing and review)

**Total**: ~3-4 weeks for complete migration

## Success Criteria

- [ ] All HE transport tests pass
- [ ] HE workflow produces same results as before restructure
- [ ] No import errors in any dependent code
- [ ] Documentation updated
- [ ] Deprecation warnings in place
- [ ] Code review approved
- [ ] CI/CD passes

---

## VS Code Task List Format

```markdown
## HE Transport Restructure Tasks

### Setup
- [ ] Create modes/he/ directory structure
- [ ] Create __init__.py files
- [ ] Create tests/ subdirectory

### Core Files Migration  
- [ ] Create solver.py with HighEnergyTransportSolver
- [ ] Migrate main_gamma_ray_loop.py → gamma_transport.py
- [ ] Migrate gamma_packet_loop.py → packet_propagation.py
- [ ] Migrate gamma_ray_interactions.py → gamma_interactions.py
- [ ] Migrate gamma_ray_grid.py → gamma_grid.py
- [ ] Migrate gamma_ray_estimators.py → estimators.py

### Packet Source
- [ ] Migrate high_energy.py → modes/he/packet_source.py
- [ ] Move or reference GXPacket appropriately
- [ ] Update imports in modes/he/__init__.py

### Tests
- [ ] Migrate HE transport tests to modes/he/tests/
- [ ] Update test imports
- [ ] Verify all tests pass

### Workflows
- [ ] Update tardis_he_workflow.py imports
- [ ] Update workflow tests
- [ ] Run end-to-end HE workflow test

### Documentation
- [ ] Create conceptual changelog entry
- [ ] Update docstrings with new import paths
- [ ] Update notebooks using HE transport

### Deprecation
- [ ] Add deprecation warnings to old locations
- [ ] Update __init__.py files
- [ ] Document migration guide

### Final Steps
- [ ] Run full test suite
- [ ] Code review
- [ ] Merge (after PR #3427)
```

## Related PRs and Issues

- PR #3427: Monte Carlo mode architecture (must merge first)
- Related to: Future HE physics enhancements
- Dependencies: None blocking, but wait for #3427

---

**Document Version**: 1.0  
**Created**: February 2026  
**Author**: Copilot Agent  
**Status**: Proposal
