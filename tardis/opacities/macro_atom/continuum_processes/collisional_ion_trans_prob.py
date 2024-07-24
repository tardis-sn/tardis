import pandas as pd

from tardis import constants as const
from tardis.plasma.properties.base import (
    ProcessingPlasmaProperty,
    TransitionProbabilitiesProperty,
)
from tardis.plasma.properties.continuum_processes.rates import (
    IndexSetterMixin,
)

H = const.h.cgs.value

class RawCollIonTransProbs(TransitionProbabilitiesProperty, IndexSetterMixin):
    """
    Attributes
    ----------
    p_coll_ion : pandas.DataFrame, dtype float
        The unnormalized transition probabilities for
        collisional ionization.
    p_coll_recomb : pandas.DataFrame, dtype float
        The unnormalized transition probabilities for
        collisional recombination.
    cool_rate_coll_ion : pandas.DataFrame, dtype float
        The collisional ionization cooling rates of the electron gas.
    """

    outputs = ("p_coll_ion", "p_coll_recomb", "cool_rate_coll_ion")
    transition_probabilities_outputs = (
        "p_coll_ion",
        "p_coll_recomb",
        "cool_rate_coll_ion",
    )
    latex_name = (
        r"p^{\textrm{coll ion}}",
        r"p^{\textrm{coll recomb}}",
        r"C^{\textrm{ion}}",
    )

    def calculate(
        self,
        coll_ion_coeff,
        coll_recomb_coeff,
        nu_i,
        photo_ion_idx,
        electron_densities,
        energy_i,
        level_number_density,
    ):
        p_coll_ion = coll_ion_coeff.multiply(energy_i, axis=0)
        p_coll_ion = p_coll_ion.multiply(electron_densities, axis=1)
        p_coll_ion = self.set_index(p_coll_ion, photo_ion_idx, reverse=False)

        coll_recomb_rate = coll_recomb_coeff.multiply(
            electron_densities, axis=1
        )  # The full rate is obtained from this by multiplying by the
        # electron density and ion number density.
        p_recomb_deactivation = coll_recomb_rate.multiply(nu_i, axis=0) * H
        p_recomb_deactivation = self.set_index(
            p_recomb_deactivation, photo_ion_idx, transition_type=-1
        )
        p_recomb_deactivation = p_recomb_deactivation.groupby(level=[0]).sum()
        index_dd = pd.MultiIndex.from_product(
            [p_recomb_deactivation.index.values, ["k"], [0]],
            names=list(photo_ion_idx.columns) + ["transition_type"],
        )
        p_recomb_deactivation = p_recomb_deactivation.set_index(index_dd)

        p_recomb_internal = coll_recomb_rate.multiply(energy_i, axis=0)
        p_recomb_internal = self.set_index(
            p_recomb_internal, photo_ion_idx, transition_type=0
        )
        p_coll_recomb = pd.concat([p_recomb_deactivation, p_recomb_internal])

        cool_rate_coll_ion = (coll_ion_coeff * electron_densities).multiply(
            nu_i * H, axis=0
        )
        level_lower_index = coll_ion_coeff.index
        cool_rate_coll_ion = (
            cool_rate_coll_ion
            * level_number_density.loc[level_lower_index].values
        )
        cool_rate_coll_ion = self.set_index(
            cool_rate_coll_ion, photo_ion_idx, reverse=False
        )
        cool_rate_coll_ion = cool_rate_coll_ion.groupby(
            level="destination_level_idx"
        ).sum()
        ion_cool_index = pd.MultiIndex.from_product(
            [["k"], cool_rate_coll_ion.index.values, [0]],
            names=list(photo_ion_idx.columns) + ["transition_type"],
        )
        cool_rate_coll_ion = cool_rate_coll_ion.set_index(ion_cool_index)
        return p_coll_ion, p_coll_recomb, cool_rate_coll_ion


class CollRecombRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    coll_recomb_coeff : pandas.DataFrame, dtype float
        The rate coefficient for collisional recombination.
        Multiply with the electron density squared and the ion number density
        to obtain the total rate.

    Notes
    -----
    The collisional recombination rate coefficient is calculated from the
    collisional ionization rate coefficient based on the requirement of detailed
    balance.
    """

    outputs = ("coll_recomb_coeff",)
    latex_name = (r"c_{\kappa\textrm{i,}}",)

    def calculate(self, phi_ik, coll_ion_coeff):
        return coll_ion_coeff.multiply(phi_ik.loc[coll_ion_coeff.index])
