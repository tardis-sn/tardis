# Monte Carlo Estimators Separation (February 11, 2026)

## Overview

Refactored the Monte Carlo radiation field estimators from a single combined class into three separate, focused classes. This improves code clarity, maintainability, and allows independent optimization of each estimator type.

## Motivation

The original `RadiationFieldMCEstimators` class contained all estimator arrays in a single jitclass:
- Cell-level radiation field properties (mean intensity, mean frequency)
- Line-specific interaction tracking (J_blue, energy deposition)
- Continuum process tracking (photoionization, heating, cooling)

This monolithic structure made it difficult to:
- Understand which estimators were used by which functions
- Optimize individual estimator types independently
- Test specific estimator functionality in isolation
- Extend or modify one estimator type without affecting others

## Architecture Changes

### New Estimator Classes

Created three separate jitclass-based estimators:

#### 1. `EstimatorsBulk` (estimators_bulk.py)
**Purpose:** Cell-level radiation field properties

**Attributes:**
- `mean_intensity_total` (J): Mean radiation field intensity per cell
- `mean_frequency` (ν̄): Mean frequency per cell

**Usage:** Updated by `update_estimators_bulk()` in `move_r_packet()`

#### 2. `EstimatorsLine` (estimators_line.py)
**Purpose:** Line-specific interaction tracking

**Attributes:**
- `mean_intensity_blue` (J_blue): Blue-shifted intensity at each line and cell
- `energy_deposition_line`: Energy deposited via line interactions

**Usage:** Updated by `update_estimators_line()` in `trace_packet()`

#### 3. `EstimatorsContinuum` (estimators_continuum.py)
**Purpose:** Continuum process tracking

**Attributes:**
- `photo_ion_estimator`: Photoionization rate estimator
- `stim_recomb_estimator`: Stimulated recombination rate estimator
- `bf_heating_estimator`: Bound-free heating estimator
- `stim_recomb_cooling_estimator`: Stimulated recombination cooling estimator
- `ff_heating_estimator`: Free-free heating estimator
- `photo_ion_estimator_statistics`: Photoionization statistics counter

**Usage:** Updated by `update_estimators_bound_free()` in `single_packet_loop()`

### Factory Pattern

Each estimator class includes:
- **Instance factory:** `init_estimators_*()` - Allocates and returns a single estimator instance
- **List factory:** `create_estimators_*_list()` - Creates multiple instances for parallel processing
- **Increment method:** Adds another estimator's arrays element-wise (for accumulating thread results)

### Global vs Local Estimators

Implemented a clear separation for parallel Monte Carlo simulations:

**Global Estimators** (suffix: `_global`):
- Single instances that accumulate final results
- Created at the start of `montecarlo_main_loop()`
- Returned as final output

**Local Estimators** (suffix: `_local`):
- One set per thread to avoid race conditions
- Each thread works with its own local estimators
- Accumulated into global estimators after parallel loop completes

Example:
```python
# Global estimators - final accumulated results
estimators_bulk_global = init_estimators_bulk(n_cells)
estimators_line_global = init_estimators_line(tau_sobolev_shape)
estimators_continuum_global = init_estimators_continuum(gamma_shape, n_cells)

# Local estimators - one per thread
estimators_bulk_list_local = create_estimators_bulk_list(n_cells, n_threads)
estimators_line_list_local = create_estimators_line_list(tau_sobolev_shape, n_threads)
estimators_continuum_list_local = create_estimators_continuum_list(gamma_shape, n_cells, n_threads)

# After parallel loop: accumulate results
for i in range(n_threads):
    estimators_bulk_global.increment(estimators_bulk_list_local[i])
    estimators_line_global.increment(estimators_line_list_local[i])
    estimators_continuum_global.increment(estimators_continuum_list_local[i])
```

## Implementation Details

### Core Files Modified

1. **montecarlo_main_loop.py**
   - Returns three separate estimators instead of one combined
   - Uses individual factory functions (`init_estimators_bulk`, `init_estimators_line`, `init_estimators_continuum`)
   - Implements `_global`/`_local` naming convention for clarity

2. **single_packet_loop.py**
   - Accepts three separate local estimators (one set per thread)
   - Passes appropriate estimators to sub-functions

3. **r_packet_transport.py**
   - `trace_packet()` accepts `estimators_line`
   - `move_r_packet()` accepts `estimators_bulk`

4. **radfield_estimator_calcs.py**
   - Renamed functions:
     - `update_line_estimators` → `update_estimators_line`
     - `update_base_estimators` → `update_estimators_bulk`
   - Updated type hints to use specific estimator classes
   - Renamed attributes to match new naming convention:
     - `j_estimator` → `mean_intensity_total`
     - `nu_bar_estimator` → `mean_frequency`
     - `j_blue_estimator` → `mean_intensity_blue`
     - `Edotlu_estimator` → `energy_deposition_line`

5. **base.py** (MonteCarloTransportSolver)
   - Receives three estimators from `montecarlo_main_loop()`
   - Stores them separately in `transport_state`

6. **montecarlo_transport_state.py**
   - Added three attributes: `estimators_bulk`, `estimators_line`, `estimators_continuum`
   - Created backward compatibility property `radfield_mc_estimators` that dynamically combines all three
   - Updated delegating properties (`j_estimator`, `nu_bar_estimator`, `j_blue_estimator`) to access individual estimators directly

### Backward Compatibility

**Legacy Support:**
- Kept `RadiationFieldMCEstimators` in `legacy_mc_estimators.py` for backward compatibility
- Created `radfield_mc_estimators` property in `montecarlo_transport_state.py` that dynamically combines the three estimators
- HDF5 serialization continues using delegating properties, maintaining data structure

**Migration Path:**
- **Preferred:** Direct access to individual estimators
  ```python
  transport_state.estimators_bulk.mean_intensity_total
  transport_state.estimators_line.mean_intensity_blue
  transport_state.estimators_continuum.photo_ion_estimator
  ```
- **Legacy:** Property access (for backward compatibility)
  ```python
  transport_state.j_estimator
  transport_state.nu_bar_estimator
  transport_state.j_blue_estimator
  ```

### Tests Fixed

Updated test files to work with new estimator structure:
- **test_packet.py**: Updated imports, fixture, and function calls
  - Changed `update_line_estimators` → `update_estimators_line`
  - Updated estimator fixture to create `EstimatorsLine` instead of `RadiationFieldMCEstimators`
  - Updated attribute names in assertions

## Terminology Changes

- **"shells"** → **"cells"**: Updated throughout to support future multi-dimensional geometries
- **Estimator attributes:** More descriptive names (e.g., `mean_intensity_total` instead of `j_estimator`)

## Benefits

1. **Finer-Grained Separation**: Bulk vs line vs continuum estimators clearly separated
2. **Clarity of Purpose**: Each file has specific, focused responsibility
3. **Maintainability**: Easier to modify one type without affecting others
4. **Performance**: Can optimize each estimator type independently
5. **Type Safety**: Clearer type hints with focused classes
6. **Testing**: Can test each estimator type independently
7. **Factory Pattern**: Clean separation between allocation and initialization
8. **Explicit Parallel Processing**: `_global`/`_local` naming makes intent clear
9. **Memory Efficiency**: List factories allocate fresh zeros instead of copying
10. **Backward Compatibility**: Legacy code continues to work

## Documentation

Created example notebook: `docs/physics/transport/montecarlo/estimators_example.ipynb`
- Demonstrates running Monte Carlo loop directly
- Shows how to access and visualize estimators
- Provides examples of bulk, line, and continuum estimator usage

## Testing

All tests pass:
- ✅ Full TARDIS regression tests (5/5)
- ✅ Monte Carlo transport tests (55 passed, 57 skipped, 6 xfailed)
- ✅ Additional transport tests (43 passed, 40 skipped, 3 xpassed)

**Regression data unchanged** - HDF5 output structure and values remain identical.

## Commands Executed

1. Initial refactoring following 13-step plan
2. Additional refinement: Renamed estimators with `_global`/`_local` suffixes for clarity
3. Used individual factory functions instead of combined `initialize_estimator_statistics` in main loop
4. Created concise documentation notebook showing Monte Carlo estimator usage

## Files Created

- `tardis/transport/montecarlo/estimators/estimators_bulk.py`
- `tardis/transport/montecarlo/estimators/estimators_line.py`
- `tardis/transport/montecarlo/estimators/estimators_continuum.py`
- `tardis/transport/montecarlo/estimators/legacy_mc_estimators.py` (renamed from `radfield_mc_estimators.py`)
- `docs/physics/transport/montecarlo/estimators_example.ipynb`

## Files Modified

- `tardis/transport/montecarlo/montecarlo_main_loop.py`
- `tardis/transport/montecarlo/single_packet_loop.py`
- `tardis/transport/montecarlo/r_packet_transport.py`
- `tardis/transport/montecarlo/estimators/radfield_estimator_calcs.py`
- `tardis/transport/montecarlo/base.py`
- `tardis/transport/montecarlo/montecarlo_transport_state.py`
- `tardis/transport/montecarlo/estimators/__init__.py`
- `tardis/transport/montecarlo/packets/tests/test_packet.py`
- `tardis/io/model/readers/continuum_radfield_properties.py`
- `tardis/tests/fixtures/conftest.py`
- `tardis/transport/montecarlo/tests/test_montecarlo.py`

## Future Work

- Consider deprecating `RadiationFieldMCEstimators` in future release
- Add deprecation warnings to legacy property access
- Migrate remaining code to use individual estimators directly
