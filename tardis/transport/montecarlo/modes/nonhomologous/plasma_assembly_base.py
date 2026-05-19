import logging

from tardis.plasma.assembly.base import (
    PlasmaSolverFactory as _BasePlasmaSolverFactory,
)

logger = logging.getLogger(__name__)


class PlasmaSolverFactory(_BasePlasmaSolverFactory):
    def assemble(
        self,
        number_densities,
        dilute_planckian_radiation_field,
        electron_densities=None,
        **kwargs,
    ):
        """
        Assemble the plasma based on the provided parameters and settings.

        Parameters
        ----------
        number_densities : dict
            Dictionary of number densities for different species.
        dilute_planckian_radiation_field : object
            The dilute Planckian radiation field object.
        electron_densities : array-like, optional
            Optional electron densities.

        Returns
        -------
        BasePlasma
            The assembled plasma object.
        """
        return super().assemble(
            number_densities,
            dilute_planckian_radiation_field,
            time_explosion=-99.0,  # BasePlasma dependency, pass dummy value
            electron_densities=electron_densities,
            **kwargs,
        )
