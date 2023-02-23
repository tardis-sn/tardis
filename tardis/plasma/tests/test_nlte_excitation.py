
import pandas as pd
import numpy as np
from numpy.testing import assert_almost_equal



from tardis.plasma.properties.nlte_rate_equation_solver import NLTERateEquationSolver
from tardis.plasma.properties.nlte_excitation_data import NLTEExcitationData

def test_main_nlte_calculation_bound_bound(
    nlte_atomic_dataset,
):
    """
    Using a simple case of nlte_ion for HI and HeII, checks if the calculate_rate_matrix generates the correct data.
    """
    simple_excitation_species = (1,0)
    copy_atomic_dataset = nlte_atomic_dataset
    copy_atomic_dataset.levels = nlte_atomic_dataset.levels[nlte_atomic_dataset.levels.index.get_level_values("level_number")<5]
    lines_filtered = nlte_atomic_dataset.lines[nlte_atomic_dataset.lines.index.get_level_values("level_number_lower")<5]
    copy_atomic_dataset.lines = lines_filtered[lines_filtered.index.get_level_values("level_number_upper")<5]
    simple_nlte_data = NLTEExcitationData(copy_atomic_dataset.lines, simple_excitation_species)
    simple_t_electrons = [0.5646738847]
    simple_j_blues = pd.DataFrame(0.5, index=copy_atomic_dataset.lines.index, columns=["0"])
    
    simple_beta_sobolev = pd.DataFrame(0.8, index=copy_atomic_dataset.lines.index, columns=["0"])

    actual_rate_matrix = NLTERateEquationSolver.main_nlte_calculation_bound_bound(
        copy_atomic_dataset.levels,
        simple_t_electrons,
        simple_j_blues,
        simple_beta_sobolev,
        simple_excitation_species,
        simple_nlte_data,
    )
    desired_rate_matrix = [[[ 1.00000000e+00],
        [ 1.00000000e+00],
        [ 1.00000000e+00],
        [ 1.00000000e+00],
        [ 1.00000000e+00]],

       [[ 3.39792096e+09],
        [-3.46604468e+10],
        [ 1.25856407e+10],
        [ 9.80445669e+08],
        [ 2.10361778e+08]],

       [[ 5.44879921e+08],
        [ 2.82383118e+10],
        [-1.31698271e+11],
        [ 5.96565953e+10],
        [ 4.66884472e+09]],

       [[ 1.89343842e+08],
        [ 3.89484475e+09],
        [ 1.06043389e+11],
        [-3.43063731e+11],
        [ 1.80741125e+11]],

       [[ 8.89125366e+07],
        [ 1.30210860e+09],
        [ 1.29641220e+10],
        [ 2.82404633e+11],
        [-1.85627187e+11]]]

    assert_almost_equal(
        desired_rate_matrix/np.array(actual_rate_matrix), np.ones_like(desired_rate_matrix), decimal=6
    )