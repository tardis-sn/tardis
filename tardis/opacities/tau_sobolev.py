import pandas as pd

from tardis.opacities.sobolev import (
    calculate_beta_sobolev as _calculate_beta_sobolev,
)
from tardis.opacities.sobolev import (
    calculate_sobolev_line_opacity as _calculate_sobolev_line_opacity,
)
from tardis.plasma.properties.base import ProcessingPlasmaProperty


class TauSobolev(ProcessingPlasmaProperty):
    """Calculate the Sobolev optical depth plasma property.

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

    def calculate(
        self,
        lines,
        level_number_density,
        time_explosion,
        stimulated_emission_factor,
    ):
        """
        Calculate Sobolev line opacity.

        Calculates the Sobolev line opacity based on the provided parameters.

        Parameters
        ----------
        lines : pandas.DataFrame
            DataFrame containing information about spectral lines.
        level_number_density : pandas.DataFrame
            DataFrame with level number densities.
        time_explosion : astropy.units.Quantity
            Time since explosion.
        stimulated_emission_factor : float
            Factor for stimulated emission.

        Returns
        -------
        pandas.DataFrame
            Calculated Sobolev line opacity values.

        Raises
        ------
        ValueError
            If any calculated tau_sobolevs are nan or inf.
        """
        return _calculate_sobolev_line_opacity(
            lines,
            level_number_density,
            time_explosion,
            stimulated_emission_factor,
        )


class BetaSobolev(ProcessingPlasmaProperty):
    """Calculate Sobolev escape probabilities.

    Attributes
    ----------
    beta_sobolev : Numpy Array, dtype float
    """

    outputs = ("beta_sobolev",)
    latex_name = (r"\beta_{\textrm{sobolev}}",)

    def calculate(self, tau_sobolevs: pd.DataFrame) -> pd.DataFrame:
        """Return escape probabilities for the supplied optical depths."""
        return _calculate_beta_sobolev(tau_sobolevs)
