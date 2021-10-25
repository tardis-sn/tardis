import numpy as np
import pandas as pd
import astropy.units as u
import pynonthermal

import tardis.constants as const

ERG2EV = u.erg.to(u.eV)
K_B = const.k_B.to("erg / K").value


class SpencerFanoSolver:
    def __init__(
        self,
        simulation,
        energy_min=1,
        energy_max=3000,
        number_of_points=2000,
        verbose=True,
    ):
        """
        Parameters
        ----------
        simulation : tardis.simulation
            tardis simulation object
        """
        plasma = simulation.plasma
        # plasma properties assigned at each radial grid cell
        self.all_ion_populations = plasma.ion_number_density
        self.all_partition_functions = plasma.partition_function
        self.all_temperatures = simulation.iterations_t_rad[0].value

        # plasma properties that are the same across the grid
        self.atomic_levels = plasma.atomic_data.levels
        self.ions = plasma.ionization_data
        self.atomic_transition_data = plasma.atomic_data.lines

        # Spencer-Fano solver parameters
        self.energy_min = energy_min
        self.energy_max = energy_max
        self.number_of_points = number_of_points
        self.verbose = verbose

    def sf_solver_grid(self):
        """Loops through all grid cells to produce a Spencer-Fano solver for each cell

        Returns
        -------
        list
            List of pynonthermal Spencer-Fano solver objects indexed by grid cell
        """
        spencer_fano_solution = []

        for grid_index, temperature in enumerate(self.all_temperatures):
            spencer_fano_solution.append(
                self.sf_solver(
                    self.all_ion_populations[grid_index],
                    self.all_partition_functions.loc[:, grid_index],
                    temperature,
                )
            )

        return spencer_fano_solution

    def simple_solver(self):
        spencer_fano_solution = []

        ions_reset = self.all_ion_populations.reset_index(drop=False)

        for grid_index, temperature in enumerate(self.all_temperatures):
            sf = pynonthermal.SpencerFanoSolver(
                emin_ev=self.energy_min,
                emax_ev=self.energy_max,
                npts=self.number_of_points,
                verbose=self.verbose,
            )

            for Z, ionstage, n_ion in zip(
                ions_reset["atomic_number"],
                ions_reset["ion_number"],
                ions_reset[grid_index],
            ):
                if n_ion > 0:
                    sf.add_ionisation(Z, ionstage + 1, n_ion)
                    sf.add_ion_ltepopexcitation(
                        Z, ionstage + 1, n_ion, temperature=temperature
                    )

            spencer_fano_solution.append(sf)

        return spencer_fano_solution

    def sf_solver(self, ion_populations, partition_functions, temperature):
        """Creates a pynonthermal Spencer-Fano solver object
        from excitations and ionizations.

        Parameters
        ----------
        ion_populations : [type]
            Ion population in the grid cell
        partition_functions : [type]
            Partition function solutions in the grid cell
        temperature : float
            Temperature of the grid cell

        Returns
        -------
        pynonthermal.SpencerFanoSolver
            Spencer-Fano equation solver object
        """
        sf = pynonthermal.SpencerFanoSolver(
            emin_ev=self.energy_min,
            emax_ev=self.energy_max,
            npts=self.number_of_points,
            verbose=self.verbose,
        )

        energy_grid_start = sf.engrid[0]

        transitions_dict = self.transitions_calculation(
            ion_populations, partition_functions, temperature, energy_grid_start
        )

        for index, n_ion in ion_populations.iteritems():
            Z, ionstage = index
            ionstage += 1

            transition = transitions_dict[index]
            transition["collstr"] = [-1] * len(transition)
            transition["forbidden"] = [False] * len(transition)
            transition["A"] = transition["A_ul"]

            for index, row in transition.iterrows():
                xs_vec = pynonthermal.excitation.get_xs_excitation_vector(
                    sf.engrid, row
                )

                levelnumberdensity = row.lower_pop
                epsilon_trans_ev = row.epsilon_trans_ev

                sf.add_excitation(
                    Z,
                    ionstage,
                    levelnumberdensity,
                    xs_vec,
                    epsilon_trans_ev,
                    transitionkey=None,
                )

            if n_ion > 0:
                sf.add_ionisation(Z, ionstage, n_ion)

        return sf

    def compute_lte_populations(
        self, ion_populations, partition_functions, temperature
    ):
        """Computes the LTE populations of the grid cell

        Parameters
        ----------
        ion_populations : [type]
            Ion population in the grid cell
        partition_functions : [type]
            Partition function solutions in the grid cell
        temperature : float
            Temperature of the grid cell

        Returns
        -------
        pandas.DataFrame
            DataFrame of the LTE populations
        """
        population_list = []

        for name, ion in self.atomic_levels.iterrows():
            ionid = (name[0], name[1])
            if ionid in self.ions:
                Z = name[0]
                ion_stage = name[1]
                level_number = name[2]
                nnion = ion_populations[(Z, ion_stage)]

                ltepartfunc = partition_functions.sum()

                # UNITS
                nnlevel = (
                    nnion
                    / ltepartfunc
                    * np.exp(-ion.energy / K_B / temperature)
                )

                poprow = (
                    Z,
                    ion_stage,
                    level_number,
                    nnlevel,
                    nnlevel,
                    nnlevel / nnion,
                )
                population_list.append(poprow)

        lte_populations = pd.DataFrame(
            population_list,
            columns=[
                "atomic_number",
                "ion_number",
                "level",
                "n_LTE",
                "n_NLTE",
                "ion_popfrac",
            ],
        )
        return lte_populations

    def transitions_calculation(
        self,
        ion_populations,
        partition_functions,
        temperature,
        energy_grid_start,
    ):
        """Calculate dictionary of transitions in the grid cell

        Parameters
        ----------
        ion_populations : [type]
            Ion population in the grid cell
        partition_functions : [type]
            Partition function solutions in the grid cell
        temperature : float
            Temperature of the grid cell
        energy_grid_start : float
            Initial energy of the energy grid in EV

        Returns
        -------
        dict
            Dictionary of all transitions in the grid cell
        """
        lte_populations = self.compute_lte_populations(
            ion_populations, partition_functions, temperature
        )

        transitions_dict = {}

        for ion_info, value in ion_populations.iteritems():
            atomic_number, ion_number = ion_info[0], ion_info[1]

            # queries dataframe
            population_current_ion = lte_populations.query(
                "atomic_number==@atomic_number & ion_number==@ion_number"
            )
            population_dict = {
                x.level: x["n_NLTE"]
                for _, x in population_current_ion.iterrows()
            }

            ion = self.atomic_transition_data.query(
                "atomic_number == @atomic_number and ion_number == @ion_number"
            )
            # Horrible syntax, probably a better way with Carsus/plasma input
            # ground_level_no_j = ion.levels.iloc[0].levelname.split("[")[0]
            # Currently set to 4 after this in Andreas code
            # top_gm_level = ion.levels[
            #    ion.levels.levelname.str.startswith(ground_level_no_j)
            # ].index.max()

            top_gm_level = 5
            transitions_dict[(atomic_number, ion_number)] = ion.query(
                "level_number_lower <= @top_gm_level", inplace=False
            ).copy()

            current_level_data = self.atomic_levels.query(
                "ion_number == @ion_number and atomic_number == @atomic_number"
            )

            # print(
            #     "with {} transitions from lower <= {}".format(
            #         len(transitions_dict[(atomic_number, ion_number)]),
            #         top_gm_level,
            #     )
            # )

            if not transitions_dict[(atomic_number, ion_number)].empty:
                # previously here: forbidden line filter
                # need to add forbidden lines through carsus?
                transition_energy = []
                lower_g = []
                upper_g = []

                for row in transitions_dict[
                    (atomic_number, ion_number)
                ].iterrows():
                    # Probably a better way to get the lower/upper levels from the row
                    level_lower = current_level_data.query(
                        "level_number == @row[0][2]"
                    )
                    level_upper = current_level_data.query(
                        "level_number == @row[0][3]"
                    )

                    lower_g.append(level_lower.g.values)
                    upper_g.append(level_upper.g.values)

                    transition_energy.append(
                        (level_upper.energy.values - level_lower.energy.values)
                        * ERG2EV
                    )

                # transition_energy = ion.transition_energy.values
                # lower_g = ion.lower_g.values
                # upper_g = ion.upper_g.values

                transitions_dict[(atomic_number, ion_number)].eval(
                    "epsilon_trans_ev = @transition_energy",
                    inplace=True,
                )

                # print(len(lower_g), len(transition_energy))
                # print(len(transitions_dict[(atomic_number, ion_number)]))

                transitions_dict[(atomic_number, ion_number)].eval(
                    "lower_g = @lower_g",
                    inplace=True,
                )
                transitions_dict[(atomic_number, ion_number)].eval(
                    "upper_g = @upper_g",
                    inplace=True,
                )

                transitions_dict[(atomic_number, ion_number)].query(
                    "epsilon_trans_ev >= @energy_grid_start", inplace=True
                )

                # print(transitions_dict[(atomic_number, ion_number)])

                # not sure what is happening here
                # note that x.name[2] is the lower level number
                transitions_dict[(atomic_number, ion_number)][
                    "lower_pop"
                ] = transitions_dict[(atomic_number, ion_number)].apply(
                    lambda x: population_dict.get(x.name[2], 0.0),
                    axis=1,
                )

        return transitions_dict
