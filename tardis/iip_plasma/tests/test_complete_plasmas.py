import pytest

from tardis.iip_plasma.standard_plasmas import LegacyPlasmaArray


@pytest.fixture
def standard_legacy_plasma(
    t_rad,
    number_densities,
    time_explosion,
    atomic_data,
    j_blues,
    link_t_rad_t_electron,
):
    return LegacyPlasmaArray(
        number_densities, atomic_data, time_explosion, t_rad
    )
