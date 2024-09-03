from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

    from tardis.io.atom_data.base import MoleculeData
    from tardis.io.atom_data.nlte_data import NLTEData


@dataclass
class TARDISData:
    atom_data: pd.DataFrame #rename
    # atomic structure
    ionization_data: pd.DataFrame
    levels: pd.DataFrame = None
    nlte_data: NLTEData = None
    # lines
    lines: pd.DataFrame = None
    photoionization_data: pd.DataFrame = None
    photo_ion_block_references: pd.DataFrame = None
    photo_ion_unique_index: pd.DataFrame = None
    two_photon_data: pd.DataFrame = None
    decay_radiation_data: pd.DataFrame = None
    synpp_refs: pd.DataFrame = None #remove?
    zeta_data: pd.DataFrame = None
    # ChiantiCollisionData | CMFGENCollisionData
    collision_data: pd.DataFrame = None
    collision_data_temperatures: pd.DataFrame = None
    yg_data: pd.DataFrame = None
    # MacroAtomData
    transition_probability_data: pd.DataFrame = None
    block_reference_data: pd.DataFrame = None
    # additional
    selected_atomic_numbers: list = None
    lines_upper2macro_reference_idx: list = None
    lines_lower2macro_reference_idx: list = None

    # versioning
    uuid1: str = None
    md5: str = None
    version: str = None

    # STARDIS
    linelist_atoms: pd.DataFrame = None
    linelist_molecules: pd.DataFrame = None
    molecule_data: MoleculeData = None
