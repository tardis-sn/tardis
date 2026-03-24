# Monte Carlo Mode Architecture Migration

**Date**: February 2026  
**PR**: #3427

## Summary

Completed migration to mode-based Monte Carlo transport architecture, separating classic (line-only) and IIP (line + continuum) transport modes into distinct, self-contained implementations.

## Changes

### Architecture
- Created `tardis.transport.montecarlo.modes.classic/` for line-only transport
- Created `tardis.transport.montecarlo.modes.iip/` for continuum-enabled transport
- Each mode has dedicated: `montecarlo_transport.py`, `packet_propagation.py`, `rad_packet_transport.py`
- Mode-specific solvers in `modes/{classic,iip}/solver.py`

### Import Path Updates
- Old: `from tardis.transport.montecarlo.montecarlo_main_loop import montecarlo_main_loop`
- New: `from tardis.transport.montecarlo.modes.classic.montecarlo_transport import montecarlo_transport as montecarlo_main_loop`

### Classic Mode (Line-Only)
- Removed all continuum process conditionals
- Returns 4-tuple: `(v_packets_energy_hist, vpacket_tracker, estimators_bulk, estimators_line)`
- Electron scattering only, no bound-free or free-free opacities

### IIP Mode (Line + Continuum)
- Continuum processes always enabled (no conditionals)
- Full relativity always enabled (no conditionals)
- Returns 3-tuple: `(estimators_bulk, estimators_line, estimators_continuum)`
- Currently does not track virtual packets in main transport loop

### Files Updated
- `tardis/transport/montecarlo/packets/virtual_packet.py`
- `tardis/transport/montecarlo/tests/test_single_packet_loop.py`
- `tardis/transport/montecarlo/tests/test_montecarlo.py`
- `benchmarks/transport_montecarlo_main_loop.py`
- `docs/tutorials/run_montecarlo_transport.ipynb`

## Migration Details

- Old parent files (`montecarlo_main_loop.py`, `single_packet_loop.py`, `r_packet_transport.py`) have been removed
- All code now uses mode-specific implementations
- Classic mode maintains full backward compatibility

## Future Work

- IIP mode virtual packet tracking
- Performance benchmarking between modes
