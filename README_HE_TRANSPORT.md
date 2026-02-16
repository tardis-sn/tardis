# HE Transport Restructure - Documentation Index

This directory contains planning documents for restructuring the High Energy (gamma-ray) transport in TARDIS to follow the new modular architecture introduced in PR #3427.

## 📋 Documents

### 1. **HE_TRANSPORT_SUMMARY.md** (Quick Start) ⭐
**Read this first!**
- 5-minute overview
- Visual before/after diagrams
- Quick reference guide
- Perfect for understanding the big picture

### 2. **HE_TRANSPORT_VSCODE_CHECKLIST.md** (Implementation) ⭐⭐⭐
**Use this for implementation!**
- Ready-to-copy VS Code task list
- 10 phases with detailed checkboxes
- Step-by-step instructions
- ~280 actionable items
- **This is what you copy into VS Code**

### 3. **HE_TRANSPORT_RESTRUCTURE_PLAN.md** (Full Details)
**Reference for deep dive**
- Complete technical plan (396 lines)
- Architecture decisions and rationale
- Open questions needing team input
- Risk analysis
- Timeline estimates
- Success criteria

## 🚀 How to Use These Documents

### For Quick Understanding
1. Read `HE_TRANSPORT_SUMMARY.md`
2. Look at the before/after diagrams
3. Understand the 10 phases

### For Implementation
1. Copy `HE_TRANSPORT_VSCODE_CHECKLIST.md` into VS Code
2. Follow phases 1-10 in order
3. Check off items as you complete them
4. Refer to full plan when you need details

### For Planning Discussions
1. Review `HE_TRANSPORT_RESTRUCTURE_PLAN.md`
2. Discuss open questions with team
3. Make architectural decisions
4. Update checklist based on decisions

## 📊 What's Being Restructured?

### Current State (Scattered)
High Energy transport code is spread across:
- `tardis/energy_input/` - Mixed physics and transport
- `tardis/energy_input/transport/` - Transport implementation
- `tardis/transport/montecarlo/packet_source/` - Packet source
- `tardis/workflows/high_energy/` - Workflow

### Target State (Consolidated)
All HE transport in one place:
```
tardis/transport/montecarlo/modes/he/
├── solver.py                    # Entry point (like classic/iip)
├── gamma_transport.py           # Main transport loop
├── packet_propagation.py        # Packet propagation
├── gamma_interactions.py        # Physics (Compton, pair creation)
├── gamma_grid.py                # Geometry and grid
├── estimators.py                # HE-specific estimators
├── packet_source.py             # GammaRayPacketSource
├── packets.py                   # GXPacket (if moved)
└── tests/                       # All HE transport tests
```

## 🎯 Goals

1. **Consistency** - HE transport follows same pattern as classic/iip modes
2. **Modularity** - Clear separation between decay physics and transport
3. **Maintainability** - All HE transport code in one logical location
4. **Testability** - Isolated, comprehensive test suite
5. **Extensibility** - Easy to add new features or modes

## ⏱️ Timeline

- **Phase 1-2**: Setup & Core (5-7 days)
- **Phase 3-4**: Packets & Tests (3-5 days)
- **Phase 5-6**: Workflows & Imports (3-4 days)
- **Phase 7-8**: Deprecation & Docs (3-4 days)
- **Phase 9-10**: Testing & Review (3-5 days)

**Total**: ~3-4 weeks

## ⚠️ Important Notes

1. **Wait for PR #3427 to merge first** - The new modes architecture must be in place
2. **This is a refactor, not new features** - Behavior should remain unchanged
3. **Maintain backward compatibility** - Use deprecation warnings, not breaking changes
4. **Test thoroughly** - HE workflow must produce identical results

## 🤔 Key Decisions Needed

These questions need team discussion (see full plan for details):

1. Should `GXPacket` move to `modes/he/` or stay in `energy_input/transport/`?
2. Which utilities stay in `energy_input/util.py` vs move to `modes/he/`?
3. Keep old files as wrappers with deprecation OR remove entirely?
4. Should HE estimators align with new estimator architecture?

## 📞 Contact

For questions about this plan:
- Review the full plan in `HE_TRANSPORT_RESTRUCTURE_PLAN.md`
- Open an issue for discussion
- Tag relevant team members in PR reviews

## 📝 Status

- **Created**: February 2026
- **Current Status**: Ready for team review
- **Next Step**: Team discussion and approval
- **Implementation**: After PR #3427 merges

---

## Quick File Access

| Document | Purpose | Size | When to Use |
|----------|---------|------|-------------|
| `HE_TRANSPORT_SUMMARY.md` | Overview | 5 KB | First read |
| `HE_TRANSPORT_VSCODE_CHECKLIST.md` | Implementation | 10 KB | **During work** |
| `HE_TRANSPORT_RESTRUCTURE_PLAN.md` | Full details | 15 KB | Reference |
| `README_HE_TRANSPORT.md` (this file) | Index | 4 KB | Navigation |

---

**Happy Restructuring! 🚀**
