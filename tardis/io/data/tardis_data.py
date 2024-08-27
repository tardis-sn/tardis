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
    levels: pd.DataFrame
    ionization_data: pd.DataFrame
    nlte_data: NLTEData
    # lines
    lines: pd.DataFrame
    photoionization_data: pd.DataFrame
    photo_ion_block_references = None
    photo_ion_unique_index = None
    two_photon_data: pd.DataFrame
    decay_radiation_data: pd.DataFrame
    synpp_refs: pd.DataFrame #remove?
    zeta_data: pd.DataFrame
    # ChiantiCollisionData | CMFGENCollisionData
    collision_data: pd.DataFrame
    collision_data_temperatures: pd.DataFrame
    yg_data: pd.DataFrame
    # MacroAtomData
    transition_probability_data: pd.DataFrame
    block_reference_data: pd.DataFrame
    # additional
    selected_atomic_numbers: list
    lines_upper2macro_reference_idx: list
    lines_lower2macro_reference_idx: list

    # versioning
    uuid1 = str
    md5 = str
    version = str

    # STARDIS
    linelist_atoms: pd.DataFrame
    linelist_molecules: pd.DataFrame
    molecule_data: MoleculeData
