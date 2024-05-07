import numpy as np
import pandas as pd
from astropy import units as u

from tardis import constants as const
from tardis.plasma.properties.base import ProcessingPlasmaProperty

SOBOLEV_COEFFICIENT = (
    (
        ((np.pi * const.e.gauss**2) / (const.m_e.cgs * const.c.cgs))
        * u.cm
        * u.s
        / u.cm**3
    )
    .to(1)
    .value
)


def calculate_sobolev_line_opacity(
    lines,
    level_number_density,
    time_explosion,
    stimulated_emission_factor,
):
    tau_sobolevs = (
        (lines.wavelength_cm * lines.f_lu).values[np.newaxis].T
        * SOBOLEV_COEFFICIENT
        * time_explosion.to(u.s).value
        * stimulated_emission_factor
        * level_number_density.reindex(lines.droplevel(-1).index).values
    )

    if np.any(np.isnan(tau_sobolevs)) or np.any(np.isinf(np.abs(tau_sobolevs))):
        raise ValueError(
            "Some tau_sobolevs are nan, inf, -inf in tau_sobolevs."
            " Something went wrong!"
        )

    return pd.DataFrame(
        tau_sobolevs,
        index=lines.index,
        columns=np.array(level_number_density.columns),
    )


class TauSobolev(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    tau_sobolev : Pandas DataFrame, dtype float
          Sobolev optical depth for each line. Indexed by line.
          Columns as zones.
    """

    outputs = ("tau_sobolevs",)
    latex_name = (r"\tau_{\textrm{sobolev}}",)
    latex_formula = (
        r"\dfrac{\pi e^{2}}{m_{e} c}f_{lu}\lambda t_{exp}\
        n_{lower} \Big(1-\dfrac{g_{lower}n_{upper}}{g_{upper}n_{lower}}\Big)",
    )

    def __init__(self, plasma_parent):
        super(TauSobolev, self).__init__(plasma_parent)
        self.sobolev_coefficient = (
            (
                ((np.pi * const.e.gauss**2) / (const.m_e.cgs * const.c.cgs))
                * u.cm
                * u.s
                / u.cm**3
            )
            .to(1)
            .value
        )

    def calculate(
        self,
        lines,
        level_number_density,
        time_explosion,
        stimulated_emission_factor,
    ):
        return calculate_sobolev_line_opacity(
            lines,
            level_number_density,
            time_explosion,
            stimulated_emission_factor,
        )
