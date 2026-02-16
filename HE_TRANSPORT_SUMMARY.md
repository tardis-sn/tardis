# HE Transport Restructure - Quick Summary

## What This Is About

PR #3427 introduces a new modular architecture for TARDIS Monte Carlo transport by creating separate "modes":
- `modes/classic/` - Standard line-only transport
- `modes/iip/` - Transport with continuum processes

**This plan extends that architecture to High Energy (gamma-ray) transport.**

## The Problem

High Energy transport code is currently scattered across:
```
❌ tardis/energy_input/main_gamma_ray_loop.py
❌ tardis/energy_input/transport/gamma_packet_loop.py  
❌ tardis/energy_input/transport/gamma_ray_interactions.py
❌ tardis/energy_input/transport/gamma_ray_grid.py
❌ tardis/transport/montecarlo/packet_source/high_energy.py
```

## The Solution

Consolidate HE transport into a new mode:
```
✅ tardis/transport/montecarlo/modes/he/
   ├── solver.py                    # Entry point
   ├── gamma_transport.py           # Main loop
   ├── packet_propagation.py        # Packet loop
   ├── gamma_interactions.py        # Physics
   ├── gamma_grid.py                # Geometry
   ├── estimators.py                # HE estimators
   ├── packet_source.py             # Packet source
   └── tests/                       # All tests
```

## Files to Use

1. **`HE_TRANSPORT_RESTRUCTURE_PLAN.md`** (396 lines)
   - Complete detailed plan with rationale
   - Architecture decisions
   - Open questions
   - Timeline estimates

2. **`HE_TRANSPORT_VSCODE_CHECKLIST.md`** (278 lines) ⭐ **USE THIS!**
   - Ready-to-use VS Code task list
   - 10 phases with checkboxes
   - Clear step-by-step instructions
   - Can be copied directly into VS Code

## Quick Start

1. **Wait for PR #3427 to merge** (new modes architecture must be in place first)
2. **Copy the checklist** from `HE_TRANSPORT_VSCODE_CHECKLIST.md` into VS Code
3. **Follow the phases** in order
4. **Check off items** as you complete them

## Key Phases

### Phase 1-2: Setup & Core Migration (5-7 days)
- Create directory structure
- Migrate main transport functions

### Phase 3-4: Packets & Tests (3-5 days)  
- Move packet source and GXPacket
- Migrate all tests

### Phase 5-6: Workflows & Imports (3-4 days)
- Update HE workflow
- Fix all import paths

### Phase 7-8: Deprecation & Docs (3-4 days)
- Add deprecation warnings
- Update documentation

### Phase 9-10: Testing & Review (3-5 days)
- Full test suite
- Code review and merge

**Total Time: ~3-4 weeks**

## Benefits

✅ **Consistency** - HE transport follows same pattern as classic/iip  
✅ **Modularity** - Clear separation of concerns  
✅ **Maintainability** - All HE code in one place  
✅ **Testability** - Isolated test suite  
✅ **Extensibility** - Easy to add new features  

## What Stays in energy_input/

Keep the physics and data, move the transport:

**KEEP** in `energy_input/`:
- `gamma_ray_channel.py` - Decay channel physics
- `decay_radiation.py` - Radioactive decay data  
- `energy_source.py` - Energy source abstractions
- `nuclear_energy_source.py` - Nuclear calculations
- `samplers.py` - Statistical samplers

**MOVE** to `modes/he/`:
- Transport loops and packet propagation
- Interactions (Compton, pair creation)
- Geometry and grid handling
- Packet source

## Visual Structure

### Before (Current)
```
tardis/
├── energy_input/                    ← Mixed: physics + transport
│   ├── main_gamma_ray_loop.py      
│   ├── gamma_ray_transport.py      
│   ├── gamma_ray_estimators.py     
│   └── transport/                  
│       ├── gamma_packet_loop.py    
│       ├── gamma_ray_interactions.py
│       └── gamma_ray_grid.py       
│
└── transport/montecarlo/
    └── packet_source/
        └── high_energy.py           ← Separated from main code
```

### After (Proposed)
```
tardis/
├── energy_input/                    ← Only physics & data
│   ├── gamma_ray_channel.py       ✅ Decay physics
│   ├── decay_radiation.py         ✅ Decay data
│   ├── energy_source.py           ✅ Abstractions
│   └── samplers.py                ✅ Samplers
│
└── transport/montecarlo/
    └── modes/
        ├── classic/                 ✅ Line-only
        ├── iip/                     ✅ Line + continuum
        └── he/                      ✅ Gamma-ray (NEW!)
            ├── solver.py           
            ├── gamma_transport.py  
            ├── packet_propagation.py
            ├── gamma_interactions.py
            ├── gamma_grid.py       
            ├── estimators.py       
            ├── packet_source.py    
            └── tests/              
```

## Next Steps

1. ✅ **Done**: Created comprehensive plan
2. ✅ **Done**: Created VS Code checklist  
3. ⏳ **Next**: Review with team
4. ⏳ **Next**: Wait for PR #3427 to merge
5. ⏳ **Next**: Begin implementation following checklist

## Questions?

See the full plan in `HE_TRANSPORT_RESTRUCTURE_PLAN.md` which includes:
- Detailed phase descriptions
- Open questions and decisions needed
- Risk analysis
- Timeline estimates
- Success criteria

---

**Created**: February 2026  
**Status**: Ready for review and implementation  
**Dependencies**: PR #3427 must merge first
