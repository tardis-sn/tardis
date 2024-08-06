from tardis.plasma.exceptions import PlasmaException
from tardis.plasma.properties.base import ProcessingPlasmaProperty


class BoundFreeOpacity(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    chi_bf : pandas.DataFrame, dtype float
        Bound-free opacity corrected for stimulated emission.
    """

    outputs = ("chi_bf",)
    latex_name = (r"\chi^{\textrm{bf}}",)

    def calculate(
        self,
        photo_ion_cross_sections,
        t_electrons,
        phi_ik,
        level_number_density,
        lte_level_number_density,
        boltzmann_factor_photo_ion,
    ):
        cross_section = photo_ion_cross_sections["x_sect"].values

        n_i = level_number_density.loc[photo_ion_cross_sections.index]
        lte_n_i = lte_level_number_density.loc[photo_ion_cross_sections.index]
        chi_bf = (n_i - lte_n_i * boltzmann_factor_photo_ion).multiply(
            cross_section, axis=0
        )

        num_neg_elements = (chi_bf < 0).sum().sum()
        if num_neg_elements:
            raise PlasmaException("Negative values in bound-free opacity.")
        return chi_bf
