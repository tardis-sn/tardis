from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

    from tardis.io.atom_data.base import MoleculeData
    from tardis.io.atom_data.nlte_data import NLTEData


@dataclass
class TARDISData:
    # Atomic symbols and masses. Used in Composition, LineInfo widget, Grotrian diagram widget
    atom_data: pd.DataFrame  # rename

    # ionization data
    # DIRECT PLASMA PROPERTY. Used for Van Regenmorter approx
    ionization_data: pd.DataFrame

    # fraction of recombination to ground vs excited states (ionization)
    # DIRECT PLASMA PROPERTY
    zeta_data: pd.DataFrame = None

    # atomic structure
    # DIRECT PLASMA PROPERTY
    levels: pd.DataFrame = None

    # lines
    # DIRECT PLASMA PROPERTY
    lines: pd.DataFrame = None

    # photo ionization
    # DIRECT PLASMA PROPERTY. Also uses levels, macro_atom_references. Also used in DiluteBlackBodyContinuumPropertiesSolver
    photoionization_data: pd.DataFrame = None
    # calculated in atom_data and also in ContinuumState
    photo_ion_block_references: pd.DataFrame = None

    # Used ONLY in DiluteBlackBodyContinuumPropertiesSolver
    photo_ion_unique_index: pd.DataFrame = None

    # DIRECT PLASMA PROPERTY
    two_photon_data: pd.DataFrame = None

    # code comparison property
    synpp_refs: pd.DataFrame = None

    # DIRECT PLASMA PROPERTY. Created using lines, collision_data, levels, collision_data_temperatures
    nlte_data: NLTEData = None

    # Chianti
    collision_data: pd.DataFrame = None
    # Chianti or CMFGEN. Used in NLTEData and YgData
    collision_data_temperatures: pd.DataFrame = None
    # CMFGEN. DIRECT PLASMA PROPERTY, Also uses collision_data, collision_data_temperatures, levels, macro_atom_references
    # Van Regenmorter uses ionization_data and lines
    yg_data: pd.DataFrame = None

    # selected atoms, used in plasma assembly only
    selected_atomic_numbers: list = None

    # MacroAtomData
    # DIRECT PLASMA PROPERTY. Also used in NonContinuumTransProbsMask, RawRadBoundBoundTransProbs, formal integral (make_source_function), all outside of plasma property
    macro_atom_data: pd.DataFrame = None
    macro_atom_references: pd.DataFrame = None

    # lines connected to levels (macro atom)
    # DIRECT PLASMA PROPERTY LevelIdxs2LineIdx as a pair
    lines_upper2macro_reference_idx: list = None
    lines_lower2macro_reference_idx: list = None

    # level connection to continuum, recalculated in PhotoIonizationData
    # Used directly in MCContinuumPropertiesSolver
    level2continuum_edge_idx: pd.Series = None

    # gamma rays
    decay_radiation_data: pd.DataFrame = None

    # versioning
    uuid1: str = None
    md5: str = None
    version: str = None

    # STARDIS (all optional)
    linelist_atoms: pd.DataFrame = None
    linelist_molecules: pd.DataFrame = None
    molecule_data: MoleculeData = None
