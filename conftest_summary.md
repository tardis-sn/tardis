# TARDIS `conftest.py` Summary

Pytest discovers `conftest.py` files upward from each test directory. The **root parent** for the whole package is `tardis/conftest.py`. Intermediate parents (e.g. `tardis/transport/montecarlo/conftest.py`) also apply to nested test dirs unless a closer conftest overrides or re-exports fixtures.

**Atom-data fixtures** (`atomic_data_fname`, `atomic_dataset`, NLTE variants, etc.) are defined in `tardis/tests/fixtures/atom_data.py` and registered globally via `from tardis.tests.fixtures.atom_data import `* in the root conftest.

---

## Root: `tardis/conftest.py`


|                   |                                                       |
| ----------------- | ----------------------------------------------------- |
| **Role**          | Package-wide pytest configuration and shared fixtures |
| **Reuses parent** | N/A (top-level)                                       |


**Hooks / CLI**

- `pytest_configure` — requires `tardisbase`, Astropy header, `ignore_generate` marker, regression-data plugin
- `pytest_addoption` — `--tardis-regression-data`, `--integration-tests`, `--generate-reference`, `--less-packets`
- `pytest_collection_modifyitems` — skips `ignore_generate` tests when generating reference data
- `pytest_sessionfinish` — closes montecarlo progress bars

**Plugins**

- `tardisbase.testing.regression_data.regression_data` (when `tardisbase` is installed)

**Fixtures defined here**


| Fixture                                     | Scope    | Notes                                         |
| ------------------------------------------- | -------- | --------------------------------------------- |
| `generate_reference`                        | session  | `--generate-reference` flag                   |
| `tardis_regression_path`                    | session  | `--tardis-regression-data` path               |
| `tardis_config_verysimple`                  | function | YAML dict                                     |
| `config_verysimple_for_simulation_one_loop` | session  | uses `config_verysimple`, `atomic_data_fname` |
| `config_verysimple_hydrogen_only`           | function | uses `config_verysimple`                      |
| `tardis_config_verysimple_nlte`             | function | NLTE YAML dict                                |
| `hdf_file_path`                             | session  | temp HDF path                                 |
| `example_model_file_dir`                    | session  | model reader test data                        |
| `example_configuration_dir`                 | session  | config test data dir                          |
| `config_verysimple`                         | session  | `Configuration` object                        |
| `config_rpacket_tracking`                   | session  | separate `Configuration` for r-packet tests   |
| `config_montecarlo_1e5_verysimple`          | function | montecarlo config                             |
| `simulation_verysimple`                     | session  | full sim, 4000 last packets                   |
| `simulation_verysimple_default`             | session  | full sim, default packets                     |
| `simulation_verysimple_vpacket_tracking`    | session  | virtual packet logging                        |
| `simulation_rpacket_tracking`               | session  | r-packet tracking enabled                     |
| `workflow_simple`                           | class    | `StandardTARDISWorkflow` run                  |
| `iip_atom_data`                             | function | IIP/ctardis atom data from regression path    |


**Also pulls in** (via star import): all fixtures from `tardis/tests/fixtures/atom_data.py`.

**Also imports from `tardis/tests`:** `monkeysession` from `tardis/tests/test_util.py` (registers it package-wide without a `tardis/tests/conftest.py`).

---

## `tardis/tests` — shared fixtures and integration tests (no conftest)

There is **no** `tardis/tests/conftest.py`. This directory is the package-level test hub: atom-data fixtures live in a submodule, two integration test modules define their own fixtures, and `test_util.py` holds util-test fixtures plus one fixture re-exported at root.

```
tardis/tests/
├── fixtures/atom_data.py      # 10 fixtures → star-imported by root conftest
├── test_util.py               # 3 fixtures (monkeysession exported at root)
├── test_tardis_full.py        # 1 class-scoped simulation fixture
└── test_tardis_full_formal_integral.py  # 3 module/class fixtures
```

### `fixtures/atom_data.py` — global atom / NLTE fixtures

Registered package-wide via `from tardis.tests.fixtures.atom_data import *` in root `conftest.py`. These are the **canonical** regression-data atom fixtures used across the suite.

| Fixture | Scope | Depends on (root) | Purpose |
|---|---|---|---|
| `atomic_data_fname` | session | `tardis_regression_path` | Path to default Kurucz H/He HDF |
| `atomic_dataset` | session | `atomic_data_fname` | Loaded `AtomData`, MD5-checked |
| `kurucz_atomic_data` | function | `atomic_dataset` | Deepcopy for mutating tests |
| `nlte_atomic_data_fname` | function | `tardis_regression_path` | Path to `TestNLTE_He_Ti.h5` |
| `nlte_atomic_dataset` | function | `nlte_atomic_data_fname` | NLTE atom data load |
| `nlte_atom_data` | function | `nlte_atomic_dataset` | Deepcopy for NLTE tests |
| `tardis_model_config_nlte_root` | function | `example_configuration_dir` | NLTE config, `nlte_solver="root"` |
| `tardis_model_config_nlte_lu` | function | `example_configuration_dir` | NLTE config, `nlte_solver="lu"` |
| `nlte_raw_model_root` | function | `tardis_model_config_nlte_root`, `nlte_atom_data` | `SimulationState` for root solver |
| `nlte_raw_model_lu` | function | `tardis_model_config_nlte_lu`, `nlte_atom_data` | `SimulationState` for LU solver |

**Consumers (examples):** `plasma/tests/test_nlte_solver.py`, `io/configuration/tests/test_config_reader.py`, `io/tests/test_atomic.py`, `plasma/tests/test_tardis_model_density_config.py`, `plasma/equilibrium/tests/test_heating_cooling_rates.py`, and most of the suite via `atomic_dataset` / `atomic_data_fname`.

**Overlap with root conftest:** root defines its own `config_verysimple`, `tardis_config_verysimple_nlte`, `iip_atom_data`, and simulation fixtures — atom_data.py does **not** duplicate those; it owns the regression atom-data chain only.

**Overlap with `plasma/equilibrium/tests/conftest.py`:** equilibrium defines separate `hydrogen_atomic_data_fname` (`cmfgen_H.h5`) and `new_chianti_atomic_dataset_si` — different files from `nlte_atomic_data_fname` / `atomic_data_fname`. Not duplicates, but parallel atom-data fixture patterns.

### `test_util.py`

| Fixture | Scope | Reuses root? | Exported globally? |
|---|---|---|---|
| `artis_data_dir` | function | No | No — local to this file |
| `artis_abundances_fname` | function | No (uses local `artis_data_dir`) | No |
| `monkeysession` | session | No | **Yes** — `from tardis.tests.test_util import monkeysession` in root conftest |

`monkeysession` is used by visualization widget tests (`test_line_info.py`, `test_shell_info.py`, `test_custom_abundance.py`) to mock UI dependencies.

**Redundancy:** `artis_data_dir` here duplicates `io/model/conftest.py` and several other test files (same path string). Only used locally for `test_convert_abundances_format`. Could use `io/model/conftest.py` if that conftest were visible to `tardis/tests/` (it is **not** — different subtree; pytest does not discover `io/model/conftest.py` for tests under `tardis/tests/`).

**Note:** `artis_abundances_fname` in `test_util.py` points at `io/model/artis/tests/data/artis_abundances.dat`; the same-named fixture in `io/model/readers/tests/test_model_reader.py` uses `example_model_file_dir / "artis_abundances.dat"` (readers test data). **Same name, different paths** — not interchangeable.

### `test_tardis_full.py`

| Fixture | Scope | Depends on | Purpose |
|---|---|---|---|
| `TestTransportSimple.simulation` | class | `atomic_data_fname`, `generate_reference`, `example_configuration_dir` | Full convergence + final run; syncs HDF regression data |

Two standalone test functions (`test_run_tardis_from_config_obj`, `test_run_tardis_simulation_callbacks_none`) inline-load config without fixtures.

**Near-duplicate of root fixtures:** class `simulation` overlaps `simulation_verysimple` / `simulation_one_loop` (same YAML, convergence + final run) but adds `RegressionData` HDF sync and class scope. Could potentially compose on root `config_verysimple` + `atomic_data_fname` instead of re-loading YAML inline.

### `test_tardis_full_formal_integral.py`

| Fixture | Scope | Depends on | Purpose |
|---|---|---|---|
| `base_config` | module (param: `downbranch` / `macroatom`) | `example_configuration_dir` | Verysimple config + formal-integral spectrum settings |
| `config` | module (param: interpolate_shells `-1` / `30`) | `base_config` | Adds shell-interpolation variant |
| `TestTransportSimpleFormalIntegral.simulation` | class | `config`, `atomic_data_fname` | Formal-integral integration run |

**Near-duplicate of root:** `base_config` starts from the same verysimple YAML as root `config_verysimple` but mutates line mode, packet counts, and spectrum method — intentionally specialized for formal-integral regression. Local `config` name shadows the generic pattern used in `simulation/tests/` and `gui/tests/` (those load plain verysimple; this one is parameterized).

### `tardis/tests` summary: redundancy vs root

| Item | Verdict |
|---|---|
| `fixtures/atom_data.py` | **Correct home** for shared atom/NLTE fixtures; properly re-exported at root |
| `monkeysession` in `test_util.py` | **Unusual** — fixture defined in a test module, imported into root conftest; would be cleaner in `fixtures/` or root conftest |
| `artis_data_dir` in `test_util.py` | **Local duplicate** of path also in `io/model/conftest.py` (not visible from this test dir) |
| `TestTransportSimple.simulation` | **Near-duplicate** of root simulation fixtures; adds regression HDF wiring |
| `base_config` / `config` in formal integral test | **Appropriately specialized** — parameterized formal-integral variants, not a plain duplicate of `config_verysimple` |

### Recommended changes for `tardis/tests`

1. Add `tardis/tests/conftest.py` (or move `monkeysession` + `artis_*` helpers to `fixtures/`) so shared helpers are not split between `test_util.py` and root imports.
2. Keep `fixtures/atom_data.py` as-is — it is the right pattern.
3. Consider having `test_tardis_full.py` build on root `config_verysimple` instead of re-loading YAML in the class fixture.
4. Rename formal-integral `config` fixture to e.g. `formal_integral_config` to avoid confusion with generic `config` fixtures elsewhere.

---

## Subdirectory conftests

### `tardis/plasma/properties/tests/conftest.py`


|                |                                                                                                                                                 |
| -------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| **Reuses root** | No explicit reuse — no fixtures defined                                                                                                         |
| **Notes**      | Import-only stub (numpy, pandas, `tardis.iip_plasma.properties`, `AtomData`). Tests rely on pytest discovering the root conftest automatically. |


---

### `tardis/plasma/equilibrium/tests/conftest.py`


|                    |                                         |
| ------------------ | --------------------------------------- |
| **Reuses root**    | **Yes** — `tardis_regression_path`      |
| **Local fixtures** | Rate-matrix / ionization solver helpers |


**Root fixtures consumed**

- `tardis_regression_path` → `new_chianti_atomic_dataset_si`, `hydrogen_atomic_data_fname`

**Fixtures defined locally**

- `mock_photoionization_cross_sections`, `new_chianti_atomic_dataset_si`, `hydrogen_atomic_data_fname`
- `radiative_rate_solver`, `collisional_rate_solver`, `rate_solver_list`
- `collisional_simulation_state`, `photoionization_rate_solver`, `collisional_ionization_rate_solver`
- `rate_matrix_solver`, `mock_boltzmann_factor`

**Also depends on** `radiative_transitions`, which is **not defined** in this conftest or the root conftest (likely missing fixture or defined elsewhere at runtime).

---

### `tardis/transport/montecarlo/conftest.py`


|                 |                                                |
| --------------- | ---------------------------------------------- |
| **Reuses root** | **Yes** — `tardis_config_verysimple_nlte`      |
| **Applies to**  | All tests under `tardis/transport/montecarlo/` |


**Fixtures defined**

- `continuum_config` — builds `Configuration` from `tardis_config_verysimple_nlte`, sets continuum species

---

### `tardis/transport/montecarlo/tests/conftest.py`


|                   |                                                                                     |
| ----------------- | ----------------------------------------------------------------------------------- |
| **Reuses root**   | **Yes** — `config_montecarlo_1e5_verysimple`, `config_verysimple`, `atomic_dataset` |
| **Also inherits** | `continuum_config` from `montecarlo/conftest.py`                                    |


**Fixtures defined**

- `montecarlo_transport_config` ← `config_montecarlo_1e5_verysimple`
- `nb_simulation_verysimple` ← `config_verysimple`, `atomic_dataset`
- `verysimple_opacity_state`, `verysimple_numba_radial_1d_geometry`, `verysimple_numba_nonhomologous_geometry`
- `verysimple_time_explosion`, `verysimple_vpacket_collection`, `verysimple_3vpacket_collection`
- `verysimple_packet_collection`, `packet`, `static_packet`
- `set_seed_fixture`, `random_call_fixture` (numba-wrapped)

---

### `tardis/transport/montecarlo/packets/tests/conftest.py`


|                   |                                                                                  |
| ----------------- | -------------------------------------------------------------------------------- |
| **Reuses parent** | **Yes — explicit star import** from `tardis.transport.montecarlo.tests.conftest` |
| **Reuses root**   | Indirectly via imported montecarlo test fixtures                                 |


No local fixtures. Re-exports everything from `montecarlo/tests/conftest.py` (which itself depends on root fixtures).

---

### `tardis/transport/montecarlo/packet_source/tests/conftest.py`


|                   |                                                                                  |
| ----------------- | -------------------------------------------------------------------------------- |
| **Reuses parent** | **Yes — explicit star import** from `tardis.transport.montecarlo.tests.conftest` |
| **Reuses root**   | Indirectly via imported montecarlo test fixtures                                 |


**Local fixtures**

- `simple_weighted_packet_source` — `BlackBodyWeightedSource` for packet-source tests

---

### `tardis/transport/montecarlo/estimators/tests/conftest.py`


|                 |                                                                                                                                                                          |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Reuses root** | **Yes** — `tardis_config_verysimple_nlte` (same as `montecarlo/conftest.py`)                                                                                             |
| **Notes**       | **Duplicate** of `tardis/transport/montecarlo/conftest.py`; tests in this dir also inherit the parent `montecarlo/conftest.py` automatically, so this file is redundant. |


**Fixtures defined**

- `continuum_config` (identical to parent montecarlo conftest)

---

### `tardis/iip_plasma/tests/conftest.py`


|                 |                                                                                                                                                       |
| --------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Reuses root** | **No** — defines its own `atomic_data` from bundled `chianti_he_db.h5` (shadows root `atomic_data` name if both were in scope; here it is local-only) |


Large local fixture graph for IIP plasma property unit tests (~40 fixtures): model inputs (`number_of_cells`, `abundance`, `density`, …), atomic properties (`levels`, `lines`, `ionization_data`, …), partition functions, ion/level populations, radiative properties (`tau_sobolev`, `beta_sobolev`, `transition_probabilities`, …).

Does **not** import or star-import from root conftest.

---

### `tardis/workflows/high_energy/tests/conftest.py`


|                 |                                                    |
| --------------- | -------------------------------------------------- |
| **Reuses root** | **No** — self-contained HE workflow config loaders |


**Fixtures defined**

- `he_test_config_dir`, `he_test_config_yaml`, `he_test_config_csvy`
- `he_test_configs` (parameterized over four merger/det config YAMLs)

---

### `tardis/visualization/conftest.py`


|                 |                                                                      |
| --------------- | -------------------------------------------------------------------- |
| **Reuses root** | **Yes** — `config_verysimple`, `atomic_dataset`, `atomic_data_fname` |


**Fixtures defined**

- `simulation_simple_tracked` ← `config_verysimple`, `atomic_dataset` (r-packet tracking, reduced packets)
- `workflow_simple_tracked` ← `config_verysimple`, `atomic_data_fname` (mirrors root `workflow_simple` pattern with tracking)

---

### `tardis/opacities/tests/conftest.py`


|                 |                                                 |
| --------------- | ----------------------------------------------- |
| **Reuses root** | **Yes** — `config_verysimple`, `atomic_dataset` |


**Fixtures defined**

- `nb_simulation_verysimple` ← same pattern as `montecarlo/tests/conftest.py` (10 iterations, no `run_final`)

---

### `tardis/spectrum/formal_integral/tests/conftest.py`


|                 |                                   |
| --------------- | --------------------------------- |
| **Reuses root** | **Yes** — `simulation_verysimple` |


**Fixtures defined**

- `simulation_verysimple_opacity_state` ← wraps root simulation in `opacity_state_initialize`

---

### `tardis/model/conftest.py`


|                 |        |
| --------------- | ------ |
| **Reuses root** | **No** |


**Fixtures defined**

- `example_csvy_file_dir` — path to `tardis/model/tests/data`

---

### `tardis/io/model/conftest.py`


|                 |        |
| --------------- | ------ |
| **Reuses root** | **No** |


**Fixtures defined**

- `artis_data_dir` — path to `tardis/io/model/artis/tests/data`

---

### `tardis/energy_input/transport/tests/conftest.py`


|                 |        |
| --------------- | ------ |
| **Reuses root** | **No** |


**Fixtures defined**

- `basic_gamma_ray` — `GXPacket` with fixed position, direction, energy, frequency

---

## `tardis/energy_input/tests` — fixtures in test modules (no conftest)

There is **no** `tardis/energy_input/tests/conftest.py`. All 19 fixtures in this directory live directly in test files. Two other test modules (`test_energy_source.py`, `test_util.py`) have no fixtures — only parametrized unit tests.

### Overview

| File | Fixtures | Reuses root? |
|---|---|---|
| `test_gamma_ray_channel.py` | 3 | Yes — `example_configuration_dir`, `atomic_dataset` |
| `test_gamma_ray_transport.py` | 3 | Yes — `example_configuration_dir`, `atomic_dataset` |
| `test_gamma_ray_packet_source_minimal.py` | 10 | Yes — `atomic_dataset` only |
| `test_energy_source.py` | 0 | — |
| `test_util.py` | 0 | — |
| `transport/tests/conftest.py` | 1 | No — `basic_gamma_ray` |

### `test_gamma_ray_channel.py`

Shared multi-isotope nebular config chain for channel/decay tests.

| Fixture | Scope | Depends on | Purpose |
|---|---|---|---|
| `gamma_ray_config` | module | `example_configuration_dir` | Loads `tardis_configv1_density_exponential_nebular_multi_isotope.yml` |
| `gamma_ray_simulation_state` | module | `gamma_ray_config`, `atomic_dataset` | Builds `SimulationState` with overridden velocity, density, time |
| `gamma_ray_test_composition` | module | `gamma_ray_simulation_state` | Returns `(isotopic_mass_fraction, cell_masses)` |

### `test_gamma_ray_transport.py`

**Duplicates the first two fixtures** from `test_gamma_ray_channel.py` verbatim, then adds a third:

| Fixture | Scope | Depends on | Purpose |
|---|---|---|---|
| `gamma_ray_config` | module | `example_configuration_dir` | **Identical** to channel test |
| `gamma_ray_simulation_state` | module | `gamma_ray_config`, `atomic_dataset` | **Identical** to channel test |
| `gamma_ray_model_state` | module | `gamma_ray_simulation_state` | Returns `(raw_isotope_abundance, cell_masses)` — similar role to `gamma_ray_test_composition` but different composition fields |

**Verdict:** `gamma_ray_config` and `gamma_ray_simulation_state` are **needless duplicates** across two sibling files. `gamma_ray_test_composition` vs `gamma_ray_model_state` are near-duplicates (both derive mass/composition data from the same simulation state).

### `test_gamma_ray_packet_source_minimal.py`

Self-contained fixture chain for `GammaRayPacketSource` regression tests. Uses the HE workflow CSVY config directly (hardcoded path to `workflows/high_energy/tests/data/tardis_he_test_config.yml`) rather than the multi-isotope YAML above.

| Fixture | Scope | Depends on | Purpose |
|---|---|---|---|
| `simulation_state` | module | `atomic_dataset` | `SimulationState.from_csvy` from HE test config |
| `total_decays_df` | module | `simulation_state` | Decay inventory via `calculate_total_decays` |
| `isotope_decay_df` | module | `total_decays_df`, `atomic_dataset` | Gamma-ray line decay table |
| `effective_time_array` | module | — | Linear time grid via `get_effective_time_array` |
| `cumulative_decays_df` | module | `simulation_state`, `atomic_dataset`, `effective_time_array` | Time-evolved cumulative decays |
| `packet_source_params` | module | above chain | Dict of kwargs for `GammaRayPacketSource` |
| `gamma_ray_source` | module | `packet_source_params` | Instantiated packet source |
| `legacy_energy_per_packet` | module | `gamma_ray_source` | Legacy energy-per-packet calculation |
| `created_packets_data_legacy` | module | `gamma_ray_source`, `cumulative_decays_df`, `legacy_energy_per_packet` | Packets via legacy energy path |
| `created_packets_data_cumulative` | module | `gamma_ray_source`, `cumulative_decays_df` | Packets via cumulative-decay energy path |

**Root reuse:** only `atomic_dataset`. Does **not** use `workflows/high_energy/tests/conftest.py` fixtures (`he_test_config_yaml`, etc.) despite loading the same YAML — **parallel setup** to the HE workflow tests.

**Verdict:** The 10-fixture chain is appropriately local to this file (deep dependency graph), but `simulation_state` overlaps conceptually with HE workflow test infrastructure and could import from `workflows/high_energy/tests/conftest.py` or a shared energy-input conftest.

### Recommended `energy_input/tests/conftest.py`

A single conftest could consolidate:

1. **Shared nebular chain** — `gamma_ray_config`, `gamma_ray_simulation_state` (used by channel + transport tests)
2. **Composition helpers** — optionally unify `gamma_ray_test_composition` / `gamma_ray_model_state` behind one fixture with a parameter, or keep separate if the differing fields are intentional
3. **Optional HE config helper** — `simulation_state` from `test_gamma_ray_packet_source_minimal.py`, or re-export from `workflows/high_energy/tests/conftest.py`

Leave `basic_gamma_ray` in `transport/tests/conftest.py` — it is transport-layer unit test data, not shared with the gamma-ray channel/transport integration tests.

---

## Quick reference: root fixture reuse


| Subdirectory / test dir | Reuses root (`tardis/conftest.py`)? | How |
|---|---|---|
| **`tardis/tests/fixtures/atom_data.py`** | **Defines** (re-exported at root) | star-imported; uses `tardis_regression_path`, `example_configuration_dir` |
| **`tardis/tests/` test modules** | Partial | integration tests use root fixtures; `test_util` fixtures mostly local |
| `plasma/properties/tests` | No (empty) | — |
| `plasma/equilibrium/tests`                 | Yes                                 | `tardis_regression_path`                                           |
| `transport/montecarlo`                     | Yes                                 | `tardis_config_verysimple_nlte`                                    |
| `transport/montecarlo/tests`               | Yes                                 | `config_`*, `atomic_dataset`                                       |
| `transport/montecarlo/packets/tests`       | Indirect                            | star-import → montecarlo/tests → root                              |
| `transport/montecarlo/packet_source/tests` | Indirect                            | star-import → montecarlo/tests → root                              |
| `transport/montecarlo/estimators/tests`    | Yes                                 | `tardis_config_verysimple_nlte` (duplicate of montecarlo conftest) |
| `iip_plasma/tests`                         | No                                  | own `atomic_data`                                                  |
| `workflows/high_energy/tests`              | No                                  | local config paths                                                 |
| `visualization`                            | Yes                                 | `config_verysimple`, `atomic_dataset`, `atomic_data_fname`         |
| `opacities/tests`                          | Yes                                 | `config_verysimple`, `atomic_dataset`                              |
| `spectrum/formal_integral/tests`           | Yes                                 | `simulation_verysimple`                                            |
| `model`                                    | No                                  | —                                                                  |
| `io/model`                                 | No                                  | —                                                                  |
| `energy_input/transport/tests`             | No                                  | —                                                                  |
| `energy_input/tests` (test modules only)   | Partial                             | `example_configuration_dir`, `atomic_dataset` in channel/transport tests |


## Test-file fixtures (beyond conftest)

Pytest also collects fixtures from `test_*.py` files. These are **module-scoped by default** (only visible to that file unless imported), but duplicate names across files still indicate copy-paste drift.

| Location | Count |
|---|---|
| `conftest.py` files (16) | 98 fixtures |
| `tardis/tests/fixtures/atom_data.py` | 10 fixtures (star-imported by root conftest) |
| `test_*.py` files with fixtures (64 files) | 205 fixtures |
| Unique names defined only in test files (not in any conftest / atom_data) | 163 |

Most test-file fixtures are **test-local** (parametrized inputs, class-scoped setup, regression CSV loaders) and are appropriately kept in the test module. The problems below are cases where a shared conftest or parent fixture already exists, or where the same name/body appears in multiple places.

---

## Redundant or needless fixture redefinitions

### High priority — exact duplicates across conftest / test files

| Fixture | Defined in | Already available from | Verdict |
|---|---|---|---|
| `nb_simulation_verysimple` | `opacities/tests/conftest.py`, `montecarlo/tests/conftest.py` | Same body in both files | **Needless duplicate** — hoist to one shared conftest (montecarlo or root) |
| `continuum_config` | `transport/montecarlo/conftest.py`, `transport/montecarlo/estimators/tests/conftest.py` | Parent montecarlo conftest is auto-discovered | **Needless duplicate** — delete estimators conftest or remove its body |
| `artis_data_dir` | `io/model/conftest.py`, `io/model/artis/tests/test_artis_readers.py`, `model/tests/test_base.py`, **`tests/test_util.py`** | `io/model/conftest.py` not visible to `tardis/tests/`; util copy is for local abundances test | **Local duplicate** in 4 places — consider `tardis/tests/conftest.py` or shared `fixtures/paths.py` |
| `artis_abundances_fname` | `tests/test_util.py`, `io/model/readers/tests/test_model_reader.py` | **Different data paths** (artis vs readers test data) | Same name, different files — rename one for clarity |
| `gamma_ray_config` | `energy_input/tests/test_gamma_ray_channel.py`, `test_gamma_ray_transport.py` | Identical body in both files | **Needless duplicate** — move to `energy_input/tests/conftest.py` |
| `gamma_ray_simulation_state` | same two energy_input test files | Identical body in both files | **Needless duplicate** — move to shared conftest |
| `gamma_ray_test_composition` / `gamma_ray_model_state` | `test_gamma_ray_channel.py`, `test_gamma_ray_transport.py` | Both derive mass/composition from same `gamma_ray_simulation_state` | **Near-duplicate** — same conftest, possibly one fixture with param |
| `simulation_state` (HE CSVY) | `energy_input/tests/test_gamma_ray_packet_source_minimal.py` | Overlaps `workflows/high_energy/tests/conftest.py` config loading | **Parallel setup** — could reuse HE workflow conftest fixtures |
| `formal_integral_geometry` | `spectrum/formal_integral/tests/test_cuda_formal_integral.py`, `test_numba_formal_integral.py` | Same parameterized geometry setup | **Likely duplicate** — could live in `formal_integral/tests/conftest.py` |
| `time_explosion` | `formal_integral/tests/test_cuda_formal_integral.py`, `test_numba_formal_integral.py` | Same constant return value | **Likely duplicate** across sibling formal-integral tests |
| `plotter`, `plotter_from_workflow`, `generate_plot_mpl_hdf`, etc. | Duplicated between `visualization/tools/tests/test_liv_plot.py` and `test_sdec_plot.py` | Same SDECPlotter wrapper pattern in both files | **Needless duplication** — shared visualization test helper conftest would help |

### Medium priority — near-duplicates (same intent, different implementation)

| Fixture / pattern | Locations | Notes |
|---|---|---|
| `config` | `simulation/tests/test_simulation.py`, `simulation/tests/test_convergence.py`, `gui/tests/test_gui.py`, `plasma/tests/test_complete_plasmas.py`, **`tests/test_tardis_full_formal_integral.py`** | All load verysimple YAML from `example_configuration_dir` or a hardcoded path. Root already exposes `config_verysimple`. Formal-integral version is parameterized (line mode + shell interpolation). |
| `simulation` / `simulation_one_loop` | **`tests/test_tardis_full.py`** (class fixture), `simulation/tests/test_simulation.py`, `gui/tests/test_gui.py` | Full test uses inline YAML load + regression HDF; overlaps root `simulation_verysimple` |
| `workflow_simple` vs `workflow_simple_tracked` | Root `conftest.py`, `visualization/conftest.py` | Nearly identical `StandardTARDISWorkflow` setup; tracked variant adds `track_rpacket=True`. Could be one fixture with a parameter or a small wrapper. |
| `simulation_simple_tracked` vs `simulation_rpacket_tracking` | `visualization/conftest.py`, root `conftest.py` | Both run tracked r-packet sims from `config_verysimple` with reduced packets. Different scopes and entry points (`run_tardis` vs `Simulation.from_config`) but overlapping purpose. |
| `time_explosion` | `iip_plasma/tests/conftest.py` (19 days), `transport/montecarlo/packets/tests/test_packet.py` (5.2e7 scalar) | Same name, unrelated values — **not duplicates**, but the packet test name shadows the IIP fixture if scopes ever overlap. Local rename in packet tests would reduce confusion. |

### Name collisions across unrelated test modules (shadowing risk)

These share a fixture name but serve different test domains. Pytest only collides if both are visible to the same test; still worth documenting:

| Name | Colliding locations |
|---|---|
| `levels`, `lines`, `ionization_data` | `iip_plasma/tests/conftest.py` vs `io/tests/test_atomic.py` |
| `artis_abundances_fname` | `io/model/readers/tests/test_model_reader.py` vs `tests/test_util.py` |

These are **scoped to separate directories** today, so they work, but the shared names make refactoring error-prone.

### Broken / inconsistent references

| Issue | Detail |
|---|---|
| `radiative_transitions` | Required by `plasma/equilibrium/tests/conftest.py` (`collisional_rate_solver`, `radiative_rate_solver`) but **never defined** anywhere in the repo |
| `new_chianti_atomic_dataset` vs `new_chianti_atomic_dataset_si` | Conftest defines `_si`; `test_radiative_rates.py` requests `new_chianti_atomic_dataset` (no `_si`) — **name mismatch**, likely a bug or stale rename |

### Appropriate test-file fixtures (not redundant)

These are correctly kept in test modules rather than conftest:

- **`plasma/equilibrium/tests/test_heating_cooling_rates.py`** — 10+ CSV regression-data loaders (`level_population`, `ion_population`, thermal estimators); only used in that file.
- **`plasma/tests/test_nlte_solver.py`** — synthetic index/rate-matrix fixtures for NLTE solver unit tests.
- **`workflows/tests/test_iip_workflow.py`** — ctardis comparison fixtures (`iip_plasma`, `iip_plasma_nlte_init`, …); workflow-specific, correctly local.
- **`workflows/high_energy/tests/test_tardis_he_workflow.py`** — session-scoped workflow run fixtures building on `he_test_configs` from conftest; appropriate layering.
- **`workflows/tests/test_workflows.py`** — workflow variants (`v_inner_workflow`, `simple_workflow_one_loop`) correctly consume root `config_verysimple_for_simulation_one_loop`.
- **`transport/montecarlo/packets/tests/test_rpacket_last_interaction_tracker.py`** — derived dataframe slices from `nb_simulation_verysimple`; test-specific extraction logic.
- **`plasma/tests/test_tardis_model_density_config.py`** — density-model config chain using root `example_model_file_dir` and `kurucz_atomic_data`; could move to conftest but not duplicated elsewhere.
- **Class-scoped fixtures** in `visualization/tools/tests/test_sdec_plot.py`, `test_liv_plot.py`, `model/tests/test_base.py` — test-class setup; conventional pytest pattern.

### Fixtures that should stay in test files (by design)

| Pattern | Examples | Why |
|---|---|---|
| Parameterized test inputs | `invalid_coefficients` in `test_radiative_rates.py`, `formal_integral_geometry` params | Tied to specific parametrized test cases |
| `autouse` module/class setup | `to_hdf_buffer` in `test_base.py`, `test_spectrum.py`, `test_base.py` (montecarlo) | Per-module side effects |
| Test-class fixtures | `TestCompletePlasmas.config` in `test_complete_plasmas.py` | Scoped to one test class |

---

## Recommended consolidations

1. **Create `tardis/transport/montecarlo/tests/conftest.py` as single owner** of `nb_simulation_verysimple`; delete copy in `opacities/tests/conftest.py` (opacities tests already inherit montecarlo fixtures if run under package, or add explicit import).
2. **Delete** `transport/montecarlo/estimators/tests/conftest.py` — parent `montecarlo/conftest.py` already provides `continuum_config`.
3. **Add `energy_input/tests/conftest.py`** with shared `gamma_ray_config`, `gamma_ray_simulation_state`, and optionally unified composition helpers; consider wiring `test_gamma_ray_packet_source_minimal.py` to `workflows/high_energy/tests/conftest.py` instead of hardcoding the HE YAML path.
4. **Remove local `artis_data_dir`** from test files; rely on `io/model/conftest.py` (and optionally move `test_util.py` copy to a shared location).
5. **Add `visualization/tools/tests/conftest.py`** for shared `plotter` / HDF plot fixtures used by SDEC and LIV tests.
6. **Fix equilibrium naming bugs**: define `radiative_transitions` (or derive from `new_chianti_atomic_dataset_si.lines`) and align `new_chianti_atomic_dataset` naming in `test_radiative_rates.py`.
7. **Replace local `config` fixtures** in simulation/convergence tests with root `config_verysimple` where equivalent.
8. **`tardis/tests`:** move `monkeysession` (and optionally path helpers) out of `test_util.py` into `tardis/tests/conftest.py` or `fixtures/`; have `test_tardis_full.py` compose on root config/simulation fixtures; rename formal-integral `config` to avoid name collision.

---

## Notable patterns (conftest-level)

1. **Star re-export** — `packets/tests` and `packet_source/tests` explicitly `import *` from `montecarlo/tests/conftest.py` rather than relying on directory hierarchy alone.
2. **Duplicated fixtures** — `nb_simulation_verysimple` exists in both `opacities/tests/conftest.py` and `montecarlo/tests/conftest.py` with the same implementation.
3. **Duplicated conftest** — `montecarlo/conftest.py` and `montecarlo/estimators/tests/conftest.py` are identical.
4. **Shadowing** — `iip_plasma/tests/conftest.py` defines a local `atomic_data` fixture separate from the regression-data atom fixtures in the root chain.
5. **Test-file sprawl** — 205 fixtures live in test modules vs 108 in conftest/atom_data; most are intentional, but ~15 names are needless duplicates of conftest or sibling-test definitions.

