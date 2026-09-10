import numpy as np
import pandas as pd

from tardis.iip_plasma.properties.radiative_properties import (
    StimulatedEmissionFactor,
)


def test_stimulated_emission_factor_nlte_lines() -> None:
    lines = pd.DataFrame(
        index=pd.MultiIndex.from_tuples(
            [(1, 0, 0, 1), (2, 0, 0, 1)],
            names=[
                "atomic_number",
                "ion_number",
                "level_number_lower",
                "level_number_upper",
            ],
        )
    )
    level_number_density = pd.DataFrame({0: [1.0, 2.0, 1.0, 2.0]})
    g = pd.Series([1.0, 1.0, 1.0, 1.0])
    metastability = pd.Series([False, False, False, False])

    factor_module = StimulatedEmissionFactor(nlte_species=[(1, 0)])

    stimulated_emission_factor = factor_module.calculate(
        g,
        level_number_density,
        np.array([0, 2]),
        np.array([1, 3]),
        metastability,
        lines,
    )

    assert stimulated_emission_factor[0, 0] == 0.0
    assert stimulated_emission_factor[1, 0] == -1.0
