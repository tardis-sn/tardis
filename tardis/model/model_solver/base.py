from abc import ABC, abstractmethod
import pandas as pd
from astropy import units as u
import radioactivedecay as rd
from radioactivedecay.utils import Z_DICT

from tardis.model.base import Composition
from tardis.model import ModelState


class ModelSolver(ABC):
    """ """

    def __init__(self, model_state: ModelState):
        self._model_state = model_state

    def evolve_composition(self, current_time):
        """Evolve the composion of isotope and non-isotope elements in the
        simulation.

        Parameters
        ----------
        current_time : int64
            Current time of the simulation.

        Returns
        -------
        Composition : class
            The evolved compostion sate of the simulation.
        """
        atomic_mass = None
        composition = self._model_state.composition

        mass = {}
        stable_atomic_numbers = (
            composition.elemental_mass_fraction.index.to_list()
        )
        for z in stable_atomic_numbers:
            mass[z] = [
                composition.elemental_mass_fraction[z]
                for i in range(composition.elemental_mass_fraction.columns.size)
            ]
        stable_isotope_mass = pd.DataFrame(mass).T

        isotope_mass = {}
        # .decay is bad
        for atomic_number, i in composition.raw_isotope_abundance.decay(
            current_time
        ).groupby(level=0):
            i = i.loc[atomic_number]
            for column in i:
                mass = {}
                shell_abundances = i[column]
                isotopic_masses = [
                    rd.Nuclide(Z_DICT[atomic_number] + str(i)).atomic_mass
                    for i in shell_abundances.index.to_numpy()
                ]
                mass[atomic_number] = (shell_abundances * isotopic_masses).sum()
                mass[atomic_number] /= shell_abundances.sum()
                mass[atomic_number] = mass[atomic_number] * u.u.to(u.g)
                if isotope_mass.get(column) is None:
                    isotope_mass[column] = {}
                isotope_mass[column][atomic_number] = mass[atomic_number]
        isotope_mass = pd.DataFrame(isotope_mass)

        atomic_mass = pd.concat([stable_isotope_mass, isotope_mass])

        new_density = self.evolve_density(current_time)

        return Composition(
            new_density,
            composition.elemental_mass_fraction,
            atomic_mass,
            atomic_mass_unit=u.g,
        )

    @abstractmethod
    def evolve_geometry(self, current_time):
        pass

    @abstractmethod
    def evolve_density(self, current_time):
        pass

    def evolve(self, current_time):
        """Evolves the current model state

        Parameters
        ----------
        current_time : np.float64
            Current time of the simulation

        Returns
        -------
        ModelState : class
            The ModelState of the simulation
        """
        new_composition = self.evolve_composition(current_time)
        new_geometry = self.evolve_geometry(current_time)

        return ModelState(
            composition=new_composition,
            geometry=new_geometry,
            time_explosion=current_time,
        )
