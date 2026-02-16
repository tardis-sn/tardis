# HE Transport Restructure - VS Code Checklist

## Overview
Restructure High Energy (gamma-ray) transport to follow the new modular architecture from PR #3427, creating a `modes/he/` directory similar to `modes/classic/` and `modes/iip/`.

---

## Phase 1: Setup Directory Structure
- [ ] Create `tardis/transport/montecarlo/modes/he/` directory
- [ ] Create `tardis/transport/montecarlo/modes/he/__init__.py`
- [ ] Create `tardis/transport/montecarlo/modes/he/tests/` directory
- [ ] Create `tardis/transport/montecarlo/modes/he/tests/__init__.py`

---

## Phase 2: Migrate Core Transport Files

### 2.1 Create Solver
- [ ] Create `modes/he/solver.py`
  - [ ] Define `HighEnergyTransportSolver` class
  - [ ] Implement `__init__()` method
  - [ ] Implement `run()` method (entry point for HE transport)
  - [ ] Follow pattern from `modes/classic/solver.py` and `modes/iip/solver.py`

### 2.2 Migrate Main Transport Loop
- [ ] Create `modes/he/gamma_transport.py`
  - [ ] Move `run_gamma_ray_loop()` from `energy_input/main_gamma_ray_loop.py`
  - [ ] Move `get_effective_time_array()` function
  - [ ] Move `calculate_electron_number_density()` function
  - [ ] Update all imports to reference new locations

### 2.3 Migrate Packet Propagation
- [ ] Create `modes/he/packet_propagation.py`
  - [ ] Move `gamma_packet_loop()` from `energy_input/transport/gamma_packet_loop.py`
  - [ ] Preserve numba decorators and optimizations
  - [ ] Update imports

### 2.4 Migrate Interactions
- [ ] Create `modes/he/gamma_interactions.py`
  - [ ] Move `compton_scatter()` from `energy_input/transport/gamma_ray_interactions.py`
  - [ ] Move `pair_creation_packet()` function
  - [ ] Move `scatter_type()` function
  - [ ] Move `get_compton_fraction_artis()` function
  - [ ] Update imports

### 2.5 Migrate Grid/Geometry
- [ ] Create `modes/he/gamma_grid.py`
  - [ ] Move `distance_trace()` from `energy_input/transport/gamma_ray_grid.py`
  - [ ] Move `move_packet()` function
  - [ ] Update imports

### 2.6 Migrate Estimators
- [ ] Create `modes/he/estimators.py`
  - [ ] Move HE estimator classes from `energy_input/gamma_ray_estimators.py`
  - [ ] Align with new estimator architecture (`EstimatorsBulk`, `EstimatorsLine`, etc.)
  - [ ] Update imports

---

## Phase 3: Migrate Packet Source

### 3.1 Move Packet Source
- [ ] Create `modes/he/packet_source.py`
  - [ ] Move `GammaRayPacketSource` from `transport/montecarlo/packet_source/high_energy.py`
  - [ ] Move `legacy_calculate_positron_fraction()` function
  - [ ] Update imports in packet source code

### 3.2 Handle GXPacket
- [ ] Decide: Move `GXPacket` to `modes/he/packets.py` OR keep in `energy_input/transport/`
- [ ] If moving: Create `modes/he/packets.py` and move `GXPacket`, `GXPacketCollection`, `GXPacketStatus`
- [ ] If keeping: Update imports to reference existing location
- [ ] Update all imports of `GXPacket` throughout codebase

### 3.3 Update modes/he/__init__.py
- [ ] Export `HighEnergyTransportSolver`
- [ ] Export key functions/classes for external use
- [ ] Add module docstring

---

## Phase 4: Migrate Tests

### 4.1 Move HE Transport Tests
- [ ] Move tests from `energy_input/transport/tests/` to `modes/he/tests/`
  - [ ] `test_gamma_ray_grid.py`
  - [ ] `test_gamma_ray_interactions.py`
  - [ ] `conftest.py` (if HE-specific)

### 4.2 Move HE Integration Tests
- [ ] Move relevant tests from `energy_input/tests/` to `modes/he/tests/`
  - [ ] `test_gamma_ray_transport.py`
  - [ ] `test_gamma_ray_packet_source_minimal.py`
  - [ ] Any other HE-specific tests

### 4.3 Update Test Imports
- [ ] Update all test imports to reference `modes/he/` instead of `energy_input/`
- [ ] Update fixtures in `conftest.py`
- [ ] Run tests to verify migration: `pytest tardis/transport/montecarlo/modes/he/tests/`

---

## Phase 5: Update Workflows

### 5.1 Update HE Workflow
- [ ] Update `workflows/high_energy/tardis_he_workflow.py`
  - [ ] Change imports from `energy_input.*` to `transport.montecarlo.modes.he.*`
  - [ ] Update to use `HighEnergyTransportSolver` interface
  - [ ] Test workflow end-to-end

### 5.2 Update Workflow Tests
- [ ] Update `workflows/high_energy/tests/test_tardis_he_workflow.py`
  - [ ] Update imports
  - [ ] Verify tests pass
  - [ ] Add new tests if needed

---

## Phase 6: Update Imports Across Codebase

### 6.1 Find All HE Import Locations
- [ ] Search codebase for imports from:
  - [ ] `from tardis.energy_input.main_gamma_ray_loop import`
  - [ ] `from tardis.energy_input.transport.gamma_packet_loop import`
  - [ ] `from tardis.energy_input.transport.gamma_ray_interactions import`
  - [ ] `from tardis.energy_input.transport.gamma_ray_grid import`
  - [ ] `from tardis.transport.montecarlo.packet_source.high_energy import`

### 6.2 Update Each Import
- [ ] Update simulation runners
- [ ] Update notebook examples
- [ ] Update documentation
- [ ] Update any other files using HE transport

---

## Phase 7: Add Deprecation Warnings

### 7.1 Update Old File Locations
- [ ] Add deprecation warning to `energy_input/main_gamma_ray_loop.py`
- [ ] Add deprecation warning to `energy_input/transport/gamma_packet_loop.py`
- [ ] Add deprecation warning to `energy_input/transport/gamma_ray_interactions.py`
- [ ] Add deprecation warning to `energy_input/transport/gamma_ray_grid.py`
- [ ] Add deprecation warning to `transport/montecarlo/packet_source/high_energy.py`

### 7.2 Keep or Remove Old Files
- [ ] Decide: Keep old files as thin wrappers with deprecation OR remove entirely
- [ ] If keeping: Create wrapper functions that import from new location
- [ ] If removing: Delete old files after verifying no external dependencies

### 7.3 Update __init__.py Files
- [ ] Update `energy_input/__init__.py` to export only decay physics
- [ ] Update `transport/montecarlo/modes/__init__.py` to include `he` mode
- [ ] Update `transport/montecarlo/packet_source/__init__.py` if needed

---

## Phase 8: Documentation

### 8.1 Create Conceptual Changelog
- [ ] Create `docs/reference/conceptual_changelog/2026_he_transport_architecture.md`
  - [ ] Document the restructuring
  - [ ] List moved files and new locations
  - [ ] Explain benefits of new structure
  - [ ] Provide migration guide for users

### 8.2 Update Docstrings
- [ ] Update docstrings in new `modes/he/` files
- [ ] Ensure all functions have proper documentation
- [ ] Update any references to old import paths

### 8.3 Update User Documentation
- [ ] Update any user guides mentioning HE transport
- [ ] Update API documentation
- [ ] Update example notebooks using HE transport

---

## Phase 9: Final Testing and Cleanup

### 9.1 Run Test Suite
- [ ] Run all HE transport tests: `pytest tardis/transport/montecarlo/modes/he/`
- [ ] Run workflow tests: `pytest tardis/workflows/high_energy/`
- [ ] Run full test suite: `pytest tardis/`
- [ ] Fix any failing tests

### 9.2 Verify No Import Errors
- [ ] Check import statements work correctly
- [ ] Verify no circular imports
- [ ] Test in fresh Python environment

### 9.3 Code Quality
- [ ] Run linter: `ruff check tardis/transport/montecarlo/modes/he/`
- [ ] Fix any linting issues
- [ ] Run formatter if needed

### 9.4 Integration Testing
- [ ] Run end-to-end HE workflow
- [ ] Verify results match pre-restructure outputs
- [ ] Test with different configurations

---

## Phase 10: Review and Merge

### 10.1 Code Review
- [ ] Create pull request
- [ ] Request reviews from team
- [ ] Address review comments
- [ ] Update based on feedback

### 10.2 CI/CD
- [ ] Verify all CI checks pass
- [ ] Check test coverage
- [ ] Verify documentation builds correctly

### 10.3 Final Steps
- [ ] Get approval from maintainers
- [ ] Squash commits if needed
- [ ] Merge to main branch (after PR #3427 merges)
- [ ] Tag release if appropriate

---

## Summary of New File Structure

```
tardis/
├── energy_input/                        # KEEP: Decay physics & data
│   ├── gamma_ray_channel.py            # Keep: decay channel handling
│   ├── decay_radiation.py              # Keep: decay data
│   ├── energy_source.py                # Keep: energy abstractions
│   ├── nuclear_energy_source.py        # Keep: nuclear calculations
│   ├── samplers.py                     # Keep: statistical samplers
│   ├── util.py                         # Keep or move utilities
│   └── transport/                      # DEPRECATE/REMOVE
│
├── transport/montecarlo/
│   ├── modes/
│   │   ├── classic/                    # Existing
│   │   ├── iip/                        # Existing
│   │   └── he/                         # NEW ⭐
│   │       ├── __init__.py
│   │       ├── solver.py               # HE solver (entry point)
│   │       ├── gamma_transport.py      # Main HE transport loop
│   │       ├── packet_propagation.py   # Gamma packet propagation
│   │       ├── gamma_interactions.py   # Physics interactions
│   │       ├── gamma_grid.py           # Geometry/grid
│   │       ├── estimators.py           # HE estimators
│   │       ├── packet_source.py        # GammaRayPacketSource
│   │       ├── packets.py              # GXPacket (optional)
│   │       └── tests/
│   │
│   └── packet_source/
│       └── high_energy.py              # DEPRECATE
│
└── workflows/
    └── high_energy/
        └── tardis_he_workflow.py       # UPDATE imports
```

---

## Success Criteria
- [ ] All HE transport tests pass
- [ ] HE workflow produces identical results to pre-restructure
- [ ] No import errors anywhere in codebase
- [ ] Documentation updated and builds correctly
- [ ] Deprecation warnings in place for old locations
- [ ] Code review approved by maintainers
- [ ] All CI/CD checks pass

---

## Notes
- **Timing**: This restructure should be done AFTER PR #3427 merges
- **Backward Compatibility**: Maintain compatibility through deprecation warnings
- **GXPacket Decision**: Needs team decision on whether to move or keep in energy_input/
- **Duration**: Estimated 3-4 weeks for complete migration and testing
