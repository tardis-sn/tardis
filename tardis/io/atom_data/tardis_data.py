from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from tardis.io.atom_data.base import MoleculeData
    from tardis.io.atom_data.nlte_data import NLTEData


@dataclass
class TARDISData:
    atom_data: pd.DataFrame
    levels: pd.DataFrame
    lines: pd.DataFrame
    ionization_data: pd.DataFrame
    nlte_data: NLTEData
    photoionization_data: pd.DataFrame
    two_photon_data = pd.DataFrame
    decay_radiation_data = pd.DataFrame
    synpp_refs = pd.DataFrame
    zeta_data = pd.DataFrame
    # ChiantiCollisionData | CMFGENCollisionData
    collision_data: pd.DataFrame
    collision_data_temperatures: pd.DataFrame
    yg_data: pd.DataFrame
    # MacroAtomData
    transition_probability_data: pd.DataFrame
    block_reference_data: pd.DataFrame
    # STARDIS
    linelist_atoms: pd.DataFrame
    linelist_molecules: pd.DataFrame
    molecule_data: MoleculeData
