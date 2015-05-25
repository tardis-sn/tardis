import pytest
import numpy as np
import os
import tardis

from tardis.plasma.properties.atomic import Levels, Lines, LinesLowerLevelIndex, LinesUpperLevelIndex, AtomicMass, IonizationData

def test_levels_property(included_he_atomic_data):

    levels_module = Levels(None)
    selected_atoms = [2]
    levels = levels_module.calculate(included_he_atomic_data, selected_atoms)
    assert np.isclose(levels.ix[2].ix[0].ix[1]['energy'], 3.17545416e-11)

def test_lines_property(included_he_atomic_data):

    lines_module = Lines(None)
    selected_atoms = [2]
    lines = lines_module.calculate(included_he_atomic_data, selected_atoms)
    assert np.isclose(lines.ix[564954]['wavelength'], 10833.307)

def test_lines_lower_level_index_property(included_he_atomic_data):
    
    selected_atoms = [2]
    levels_module = Levels(None)
    lines_module = Lines(None)
    levels = levels_module.calculate(included_he_atomic_data, selected_atoms)
    lines = lines_module.calculate(included_he_atomic_data, selected_atoms)
    lines_lower_module = LinesLowerLevelIndex(None)
    lines_lower_level_index = lines_lower_module.calculate(levels, lines)
    assert lines_lower_level_index[9]==0

def test_lines_upper_level_index_property(included_he_atomic_data):
    
    selected_atoms = [2]
    levels_module = Levels(None)
    lines_module = Lines(None)
    levels = levels_module.calculate(included_he_atomic_data, selected_atoms)
    lines = lines_module.calculate(included_he_atomic_data, selected_atoms)
    lines_upper_module = LinesUpperLevelIndex(None)
    lines_upper_level_index = lines_upper_module.calculate(levels, lines)
    assert lines_upper_level_index[9]==30

def test_ionization_data_property(included_he_atomic_data):
    
    ionization_data_module = IonizationData(None)
    selected_atoms = [2]
    ionization_data = ionization_data_module.calculate(included_he_atomic_data, selected_atoms)
    assert np.isclose(float(ionization_data.ix[2].ix[1]),  3.9393336e-11)
