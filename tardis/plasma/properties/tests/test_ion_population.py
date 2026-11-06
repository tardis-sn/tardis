import numpy as np
import pandas as pd
from astropy import units as u

from tardis.plasma.properties.general import (
    BetaElectron,
    BetaRadiation,
    GElectron,
)
from tardis.plasma.properties.ion_population import (
    PhiSahaNebular,
    RadiationFieldCorrection,
)


def test_nebular_ionization_factor_finite_when_electrons_hotter_than_field():
    """Keep the nebular ionization factor finite for very high energies.

    Claim: When T_e exceeds T_rad, the nebular ionization factor of an ion
    with a very high ionization energy is finite, not NaN.
    Regime: One dilute shell (T_rad = 8000 K, W = 0.03, T_e = 1.1 T_rad) and
    a 10 keV ion, as for the bare-nucleus stages of Ni, whose radiation-field
    correction alone exceeds the largest double while its LTE Saha factor
    underflows to zero.
    Verification: Every ionization factor must be finite.
    """
    index = pd.MultiIndex.from_tuples(
        [(2, 0), (2, 1), (2, 2)], names=["atomic_number", "ion_number"]
    )
    partition_function = pd.DataFrame({0: [1.0, 2.0, 1.0]}, index=index)
    upper_ions = index[1:]
    ionization_data = pd.Series(
        (np.array([24.6, 1e4]) * u.eV).to_value(u.erg), index=upper_ions
    )
    zeta_data = pd.DataFrame(
        {2000.0: [0.4, 0.4], 40000.0: [0.4, 0.4]}, index=upper_ions
    )
    t_rad = np.array([8000.0])
    t_electrons = 1.1 * t_rad
    w = np.array([0.03])
    beta_rad = BetaRadiation(None).calculate(t_rad)
    beta_electron = BetaElectron(None).calculate(t_electrons)

    delta = RadiationFieldCorrection().calculate(
        w, ionization_data, beta_rad, t_electrons, t_rad, beta_electron
    )
    phi = PhiSahaNebular.calculate(
        t_rad,
        w,
        zeta_data,
        t_electrons,
        delta,
        GElectron(None).calculate(beta_rad),
        beta_rad,
        partition_function,
        ionization_data,
    )

    assert np.isfinite(phi.to_numpy()).all(), "nebular ionization factor is NaN"
