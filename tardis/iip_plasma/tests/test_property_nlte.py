import numpy as np
import pytest

from tardis.iip_plasma.standard_plasmas import LegacyPlasmaArray

pytestmark = pytest.mark.skip("Skipping tests due to old format")


def test_he_nlte_plasma(
    number_density, atomic_data, time_explosion, t_rad, w, j_blues
):
    he_nlte_plasma = LegacyPlasmaArray(
        number_densities=number_density,
        atomic_data=atomic_data,
        time_explosion=time_explosion,
        helium_treatment="recomb-nlte",
    )
    he_nlte_plasma.update_radiationfield(t_rad, w, j_blues, nlte_config=None)
    assert np.allclose(
        he_nlte_plasma.get_value("ion_number_density").loc[2].loc[1],
        number_density.loc[2],
    )
    assert (
        np.all(
            he_nlte_plasma.get_value("level_number_density")
            .loc[2]
            .loc[0]
            .loc[0]
        )
        == 0.0
    )
    assert np.allclose(
        he_nlte_plasma.get_value("level_number_density").loc[2].sum(),
        number_density.loc[2],
    )
