import pytest
from tardis.plasma.standard_plasmas import LTEPlasma

@pytest.fixture
def standard_lte_plasma_he_db(t_rad, abundance, density, time_explosion,
                              atomic_data, j_blues,
                              link_t_rad_t_electron):
    return LTEPlasma(t_rad, abundance, density, time_explosion,
                     atomic_data, j_blues, link_t_rad_t_electron)
