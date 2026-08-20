from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import pandas as pd

FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int64]
BoolArray = npt.NDArray[np.bool_]


@dataclass(frozen=True)
class ContinuumRateCoefficients:
    """Level-resolved continuum rate coefficients."""

    photoionization: FloatArray
    collisional_ionization: FloatArray
    spontaneous_recombination: FloatArray
    stimulated_recombination: FloatArray
    collisional_recombination: FloatArray


@dataclass(frozen=True)
class ContinuumCoefficientState:
    """Candidate-temperature continuum coefficients for all active levels."""

    photoionization: pd.DataFrame
    stimulated_recombination: pd.DataFrame
    spontaneous_recombination: pd.DataFrame
    collisional_ionization: pd.DataFrame


@dataclass(frozen=True)
class LevelEquationRates:
    """Density-specific rates used by one reduced level residual."""

    ionization: FloatArray
    recombination: FloatArray
    ionization_loss_matrix: FloatArray


@dataclass(frozen=True)
class BoundBoundMatrixRates:
    """Array inputs for one shell's bound-bound rate matrix."""

    number_of_levels: int
    source_level_idx: IntArray
    destination_level_idx: IntArray
    radiative_rate_coefficient: FloatArray
    collisional_rate: FloatArray
    beta_line_idx: IntArray


@dataclass(frozen=True)
class NumberDensityPerShell:
    """Absolute population information fixed for one shell."""

    hydrogen_number_density: float
    level_number_density: FloatArray
    species_level_positions: IntArray


@dataclass(frozen=True)
class SobolevInputs:
    """Line geometry required to calculate Sobolev tau and beta."""

    lines_lower_level_index: IntArray
    lines_upper_level_index: IntArray
    g_lower: FloatArray
    g_upper: FloatArray
    metastable_upper: BoolArray
    nlte_lines_mask: BoolArray
    tau_coefficient: FloatArray
    line_indices: IntArray
    line_index: pd.Index
