import logging

import numpy as np
import pandas as pd
from astropy import units as u
from scipy.interpolate import interp1d
from scipy.optimize import least_squares as lsq
from scipy.sparse import block_diag

from tardis import constants as const
from tardis.iip_plasma.continuum.base_continuum import BaseContinuum
from tardis.iip_plasma.continuum.base_continuum_data import ContinuumData
from tardis.iip_plasma.continuum.input_data import ContinuumInputData
from tardis.iip_plasma.continuum.radiative_processes import (
    RadiativeIonization,
    RadiativeRecombination,
)
from tardis.iip_plasma.standard_plasmas import LegacyPlasmaArray
from tardis.io.atom_data.parse_atom_data import parse_atom_data
from tardis.model import SimulationState
from tardis.opacities.macro_atom.macroatom_solver import (
    ContinuumMacroAtomSolver,
)
from tardis.opacities.opacity_solver import OpacitySolver
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.simulation.convergence import ConvergenceSolver
from tardis.spectrum.base import SpectrumSolver
from tardis.spectrum.formal_integral.formal_integral_solver import (
    FormalIntegralSolver,
)
from tardis.spectrum.luminosity import (
    calculate_filtered_luminosity,
)
from tardis.transport.montecarlo.base import MonteCarloTransportSolver
from tardis.transport.montecarlo.configuration import montecarlo_globals
from tardis.transport.montecarlo.progress_bars import initialize_iterations_pbar
from tardis.util.environment import Environment
from tardis.workflows.workflow_logging import WorkflowLogging

# logging support
logger = logging.getLogger(__name__)


class TypeIIPWorkflow(WorkflowLogging):
    show_progress_bars = Environment.allows_widget_display()
    enable_virtual_packet_logging = False
    log_level = None
    specific_log_level = None

    def __init__(self, configuration, csvy=False):
        """A TARDIS workflow for simulating Type IIP supernovae.

        Parameters
        ----------
        configuration : Configuration
            Configuration object for the simulation
        csvy : bool, optional
            Set true if the configuration uses CSVY, by default False
        """
        super().__init__(configuration, self.log_level, self.specific_log_level)
        self.atom_data = parse_atom_data(configuration)

        # set up states and solvers
        if csvy:
            self.simulation_state = SimulationState.from_csvy(configuration)
            assert np.isclose(
                self.simulation_state.v_inner_boundary.to(u.km / u.s).value,
                self.simulation_state.geometry.v_inner[0].to(u.km / u.s).value,
            ), (
                "If using csvy density input in the workflow, the initial v_inner_boundary must start at the first shell, see issue #3129."
            )

        else:
            self.simulation_state = SimulationState.from_config(
                configuration,
                atom_data=self.atom_data,
            )

        configuration.plasma.nlte.species = [
            (1, 0)
        ]  # Hack to force config necessary for Christian's plasma
        self.configuration = configuration

        montecarlo_globals.CONTINUUM_PROCESSES_ENABLED = True

        self.atom_data.prepare_atom_data([1], "macroatom", [(1, 0)], [(1, 0)])
        self.atom_data.continuum_data = ContinuumData(
            self.atom_data, selected_continuum_species=[(1, 0)]
        )
        # hacky thing from CTARDIS
        self.atom_data.continuum_data.photoionization_data.loc[
            (1, 0, 0), "x_sect"
        ] *= 0.0
        self.atom_data.yg_data.columns = list(
            self.atom_data.collision_data_temperatures
        )
        self.atom_data.nlte_data._init_indices()
        self.atom_data.has_collision_data = False

        elemental_number_density = (
            self.simulation_state.calculate_elemental_number_density(
                self.atom_data.atom_data.mass
            )
        )

        self.plasma_solver = LegacyPlasmaArray(
            elemental_number_density,
            self.atom_data,
            self.configuration.supernova.time_explosion.to("s"),
            nlte_config=self.configuration.plasma.nlte,
            delta_treatment=None,
            ionization_mode="nlte",  # configuration.plasma.ionization, nlte is currently not an allowed value - should be included
            excitation_mode=self.configuration.plasma.excitation,
            line_interaction_type=self.configuration.plasma.line_interaction_type,
            link_t_rad_t_electron=self.configuration.plasma.link_t_rad_t_electron
            * np.ones(self.simulation_state.geometry.no_of_shells_active),
            helium_treatment="none",
            heating_rate_data_file=None,
            v_inner=None,
            v_outer=None,
            continuum_treatment=True,
        )

        t_radiative, dilution_factor = self.initialize_radiation_field(
            self.simulation_state.geometry,
            elemental_number_density,
            self.configuration.plasma.initial_t_inner,
            self.simulation_state.dilution_factor,
        )

        radiation_field = DilutePlanckianRadiationField(
            temperature=t_radiative,
            dilution_factor=dilution_factor,
        )

        self.simulation_state.radiation_field_state = radiation_field

        j_blues = radiation_field.calculate_mean_intensity(
            self.plasma_solver.atomic_data.lines.nu.values
        )

        # Initialize NLTE

        self.plasma_solver.update_radiationfield(
            t_radiative,
            dilution_factor,
            pd.DataFrame(
                j_blues, index=self.plasma_solver.atomic_data.lines.index
            ),
            self.configuration.plasma.nlte,
            initialize_nlte=True,
            n_e_convergence_threshold=0.05,
        )

        self.base_continuum = BaseContinuum(
            plasma_array=self.plasma_solver,
            atom_data=self.atom_data,
            ws=self.simulation_state.dilution_factor,
            radiative_transition_probabilities=self.plasma_solver.transition_probabilities,
            estimators=None,
        )

        # After initializing NLTE
        line_interaction_type = configuration.plasma.line_interaction_type
        continuum_interactions = configuration.plasma.continuum_interaction
        if "H I" not in continuum_interactions.species:
            raise ValueError(
                "Continuum interactions for 'H I' must be included for the IIP workflow. Check plasma.continuum_interaction.species in the configuration."
            )

        self.opacity_solver = OpacitySolver(
            line_interaction_type,
            configuration.plasma.disable_line_scattering,
        )

        self.macro_atom_solver = ContinuumMacroAtomSolver(
            self.atom_data.levels,
            self.atom_data.lines,
            self.atom_data.photoionization_data,
            line_interaction_type=line_interaction_type,
        )

        self.transport_state = None
        self.transport_solver = MonteCarloTransportSolver.from_config(
            configuration,
            packet_source=self.simulation_state.packet_source,
            enable_virtual_packet_logging=self.enable_virtual_packet_logging,
        )

        # Luminosity filter frequencies
        self.luminosity_nu_start = (
            configuration.supernova.luminosity_wavelength_end.to(
                u.Hz, u.spectral()
            )
        )

        if u.isclose(
            configuration.supernova.luminosity_wavelength_start, 0 * u.angstrom
        ):
            self.luminosity_nu_end = np.inf * u.Hz
        else:
            self.luminosity_nu_end = (
                const.c / configuration.supernova.luminosity_wavelength_start
            ).to(u.Hz)

        # montecarlo settings
        self.total_iterations = int(configuration.montecarlo.iterations)

        self.real_packet_count = int(configuration.montecarlo.no_of_packets)

        final_iteration_packet_count = (
            configuration.montecarlo.last_no_of_packets
        )

        if (
            final_iteration_packet_count is None
            or final_iteration_packet_count < 0
        ):
            final_iteration_packet_count = self.real_packet_count

        self.final_iteration_packet_count = int(final_iteration_packet_count)

        self.virtual_packet_count = int(
            configuration.montecarlo.no_of_virtual_packets
        )

        # spectrum settings
        self.integrated_spectrum_settings = configuration.spectrum.integrated
        self.spectrum_solver = SpectrumSolver.from_config(configuration)

        # Convergence settings
        self.consecutive_converges_count = 0
        self.converged = False
        self.completed_iterations = 0
        self.luminosity_requested = (
            configuration.supernova.luminosity_requested.cgs
        )

        # Convergence solvers
        self.convergence_strategy = (
            configuration.montecarlo.convergence_strategy
        )

        self.convergence_solvers = {}
        self.convergence_solvers["t_radiative"] = ConvergenceSolver(
            self.convergence_strategy.t_rad
        )
        self.convergence_solvers["dilution_factor"] = ConvergenceSolver(
            self.convergence_strategy.w
        )
        self.convergence_solvers["t_inner"] = ConvergenceSolver(
            self.convergence_strategy.t_inner
        )

    @staticmethod
    def initialize_radiation_field(
        geometry_state, number_density, initial_t_inner, ws
    ):
        r_inner = geometry_state.r_inner_active.value
        r_outer = geometry_state.r_outer_active.value
        r_middle = (
            geometry_state.r_middle_active.value
        )  # (r_outer + r_inner) / 2.0
        delta_r = r_outer - r_inner
        t_inner = initial_t_inner

        v_inner = geometry_state.v_inner_active.value
        doppler_factor = 1.0 - v_inner / const.c.cgs.value

        sigma_T = const.sigma_T.cgs.value

        N_H = number_density.loc[1].values

        # alternative tau calculation from ctardis
        # tau_e_shell = sigma_T * delta_r * N_H

        # tau = tau_e_shell - (1 - doppler_factor) * tau_e_shell ** (
        #    2.25 * doppler_factor
        # )

        tau = sigma_T * delta_r * N_H
        tau = tau[::-1].cumsum()[::-1]
        T_eff4 = t_inner**4 / (tau[0] + 2.0 / 3.0)
        tau_middle = interp1d(r_inner, tau, fill_value="extrapolate")(r_middle)
        t_rads = (T_eff4 * (tau_middle + 2.0 / 3.0)) ** 0.25

        flat_T_start = np.where(tau_middle < 12.0)[0][0]
        t_rads[tau_middle < 12.0] = t_rads[flat_T_start]

        # Setup ws
        tau_flat = 10
        tau_geom = 0.08

        geometric_mask = tau_middle < tau_geom
        flat_mask = tau_middle > tau_flat

        a = (1.0 - tau_geom) / (np.log(tau_flat) - np.log(tau_geom))
        b = 1.0 - np.log(tau_flat) * a
        ws[flat_mask] = 1.0
        lin_mask = np.logical_and(
            np.logical_not(flat_mask), np.logical_not(geometric_mask)
        )
        ws[lin_mask] = a * np.log(tau_middle[lin_mask]) + b

        return t_rads, ws

    def get_convergence_estimates(self):
        """Compute convergence estimates from the transport state

        Returns
        -------
        dict
            Convergence estimates
        EstimatedRadiationFieldProperties
            Dilute radiation file and j_blues dataclass
        """
        estimated_radfield_properties = (
            self.transport_solver.radfield_prop_solver.solve(
                self.transport_state.estimators_bulk,
                self.transport_state.estimators_line,
                self.transport_state.time_explosion,
                self.transport_state.time_of_simulation,
                self.transport_state.geometry_state.volume,
                self.transport_state.opacity_state.line_list_nu,
            )
        )

        estimated_t_radiative = estimated_radfield_properties.dilute_blackbody_radiationfield_state.temperature
        estimated_dilution_factor = estimated_radfield_properties.dilute_blackbody_radiationfield_state.dilution_factor

        emitted_luminosity = calculate_filtered_luminosity(
            self.transport_state.emitted_packet_nu,
            self.transport_state.emitted_packet_luminosity,
            self.luminosity_nu_start,
            self.luminosity_nu_end,
        )
        absorbed_luminosity = calculate_filtered_luminosity(
            self.transport_state.reabsorbed_packet_nu,
            self.transport_state.reabsorbed_packet_luminosity,
            self.luminosity_nu_start,
            self.luminosity_nu_end,
        )

        luminosity_ratios = (
            (emitted_luminosity / self.luminosity_requested).to(1).value
        )

        estimated_t_inner = (
            self.simulation_state.t_inner
            * luminosity_ratios
            ** self.convergence_strategy.t_inner_update_exponent
        )

        logger.info(
            f"\n\tLuminosity emitted   = {emitted_luminosity:.3e}\n"
            f"\tLuminosity absorbed  = {absorbed_luminosity:.3e}\n"
            f"\tLuminosity requested = {self.luminosity_requested:.3e}\n"
        )

        self.log_plasma_state(
            self.simulation_state.t_radiative,
            self.simulation_state.dilution_factor,
            self.simulation_state.t_inner,
            estimated_t_radiative,
            estimated_dilution_factor,
            estimated_t_inner,
        )
        # ctardis does not update t_inner
        return {
            "t_radiative": estimated_t_radiative,
            "dilution_factor": estimated_dilution_factor,
            "t_inner": self.simulation_state.t_inner,
        }, estimated_radfield_properties

    def check_convergence(
        self,
        estimated_values,
    ):
        """Check convergence status for a dict of estimated values

        Parameters
        ----------
        estimated_values : dict
            Estimates to check convergence

        Returns
        -------
        bool
            If convergence has occurred
        """
        convergence_statuses = []

        for key, solver in self.convergence_solvers.items():
            current_value = getattr(self.simulation_state, key)
            estimated_value = estimated_values[key]
            no_of_shells = (
                self.simulation_state.no_of_shells if key != "t_inner" else 1
            )
            convergence_statuses.append(
                solver.get_convergence_status(
                    current_value, estimated_value, no_of_shells
                )
            )

        if np.all(convergence_statuses):
            hold_iterations = self.convergence_strategy.hold_iterations
            self.consecutive_converges_count += 1
            logger.info(
                f"Iteration converged {self.consecutive_converges_count:d}/{(hold_iterations + 1):d} consecutive "
                f"times."
            )
            return self.consecutive_converges_count >= hold_iterations + 1

        self.consecutive_converges_count = 0
        return False

    def solve_simulation_state(
        self,
        estimated_values,
    ):
        """Update the simulation state with new inputs computed from previous
        iteration estimates.

        Parameters
        ----------
        estimated_values : dict
            Estimated from the previous iterations
        """
        next_values = {}

        for key, solver in self.convergence_solvers.items():
            if (
                key == "t_inner"
                and (self.completed_iterations + 1)
                % self.convergence_strategy.lock_t_inner_cycles
                != 0
            ):
                next_values[key] = getattr(self.simulation_state, key)
            else:
                next_values[key] = solver.converge(
                    getattr(self.simulation_state, key), estimated_values[key]
                )

        self.simulation_state.t_radiative = next_values["t_radiative"]
        self.simulation_state.dilution_factor = next_values["dilution_factor"]
        self.simulation_state.blackbody_packet_source.temperature = next_values[
            "t_inner"
        ]

        return next_values

    def update_continuum_estimators(self):
        continuum_estimators = {}

        continuum_estimators["photo_ion_estimator"] = (
            self.transport_state.estimators_continuum.photo_ion_estimator
        )
        continuum_estimators["photo_ion_statistics"] = (
            self.transport_state.estimators_continuum.photo_ion_estimator_statistics
        )
        continuum_estimators["stim_recomb_estimator"] = (
            self.transport_state.estimators_continuum.stim_recomb_estimator
        )
        continuum_estimators["bf_heating_estimator"] = (
            self.transport_state.estimators_continuum.bf_heating_estimator
        )
        continuum_estimators["stim_recomb_cooling_estimator"] = (
            self.transport_state.estimators_continuum.stim_recomb_cooling_estimator
        )
        continuum_estimators["coll_deexc_heating_estimator"] = None
        continuum_estimators["ff_heating_estimator"] = (
            self.transport_state.estimators_continuum.ff_heating_estimator
        )

        j_blues_df = pd.DataFrame(
            self.transport_state.estimators_line.mean_intensity_blue,
            index=self.plasma_solver.lines.index,
        )

        continuum_estimators, j_blues_df = self.normalize_continuum_estimators(
            continuum_estimators,
            j_blues_df,
            self.transport_state.estimators_bulk.mean_intensity_total,
        )

        return continuum_estimators, j_blues_df

    def solve_plasma(self, continuum_estimators, j_blues_df):
        """Update the plasma solution with the new radiation field estimates

        Parameters
        ----------
        estimated_radfield_properties : EstimatedRadiationFieldProperties
            The radiation field properties to use for updating the plasma

        Raises
        ------
        ValueError
            If the plasma solver radiative rates type is unknown
        """
        self.plasma_solver.update_radiationfield(
            self.simulation_state.t_radiative.value,
            self.simulation_state.dilution_factor,
            j_blues_df,
            self.configuration.plasma.nlte,
            initialize_nlte=False,
            n_e_convergence_threshold=0.05,
            **continuum_estimators,
        )

    def thermal_balance_iteration(self, initial, n_e_max, nfev):
        nfev += 1
        n_e_frac = initial[::2]
        link_t_rad_t_electron = initial[1::2]

        print(f"Nfev: {nfev} \n")
        print("link:", link_t_rad_t_electron)

        pl = self.plasma_solver

        electron_densities = n_e_max * n_e_frac

        self.plasma_solver.update(
            previous_ion_number_density=pl.ion_number_density.copy(),
            previous_electron_densities=electron_densities,
            previous_beta_sobolev=pl.beta_sobolev.copy(),
            link_t_rad_t_electron=link_t_rad_t_electron,
            previous_b=pl.b,
            previous_t_electrons=pl.t_rad * link_t_rad_t_electron,
        )

        output = np.zeros(2 * len(self.plasma_solver.fractional_heating))
        frac_e_change = (
            pl.electron_densities - electron_densities
        ) / electron_densities
        n_e_frac_new = 1 - pl.electron_densities / n_e_max
        n_e_frac_change = (n_e_frac_new - (1.0 - n_e_frac)) / (1.0 - n_e_frac)

        if (
            np.logical_not(
                np.isfinite(self.plasma_solver.fractional_heating)
            ).sum()
            > 0
        ):
            print("Heating not finite\n")
        if np.logical_not(np.isfinite(frac_e_change)).sum() > 0:
            print("frac e change not finite\n")

        output[::2] = frac_e_change
        output[1::2] = self.plasma_solver.fractional_heating
        print("Frac e change:", frac_e_change)
        print("n_e_frac_change", n_e_frac_change)
        print("Heating:", self.plasma_solver.fractional_heating)
        return output

    def solve_thermal_balance(self):
        link_t_rad_t_electron_start = self.plasma_solver.link_t_rad_t_electron
        if np.array_equal(
            link_t_rad_t_electron_start,
            np.ones_like(link_t_rad_t_electron_start),
        ):
            link_t_rad_t_electron_start = (
                self.simulation_state.radiation_field_state.dilution_factor
                ** 0.25
            )
            print("Setting initial guess for link from ws:")
            print(link_t_rad_t_electron_start)

        n_e_max = (
            (
                self.plasma_solver.number_density.multiply(
                    self.plasma_solver.number_density.index.values, axis=0
                )
            )
            .sum()
            .values
        )
        n_e_frac_start = (
            self.plasma_solver.electron_densities / n_e_max
        ).values

        print("n_e_frac:", n_e_frac_start)

        initial = np.zeros(2 * len(link_t_rad_t_electron_start))
        initial[::2] = n_e_frac_start
        initial[1::2] = link_t_rad_t_electron_start
        no_shells = self.simulation_state.geometry.no_of_shells_active

        nfev = 0

        jac_sparsity = block_diag([np.ones((2, 2))] * no_shells)
        t_floor = 1500.0 * u.K
        link_floor = t_floor / self.simulation_state.t_radiative.min()
        print("Floor Link:", link_floor)

        lbound = [0.0, link_floor] * no_shells
        ubound = [1.0, 1.5] * no_shells
        self.plasma_solver.plasma_converged = False
        thermal_lsq_result = lsq(
            self.thermal_balance_iteration,
            initial,
            bounds=(lbound, ubound),
            jac_sparsity=jac_sparsity,
            xtol=1e-14,
            ftol=1e-12,
            x_scale="jac",
            verbose=1,
            max_nfev=100,
            method="trf",
            gtol=1e-14,
            args=(
                n_e_max,
                nfev,
            ),
        )
        self.plasma_solver.plasma_converged = True
        # final thermal_balance_iteration to set values in plasma
        self.thermal_balance_iteration(thermal_lsq_result.x, n_e_max, nfev)

        ion_ratio = (
            self.plasma_solver.ion_number_density.loc[(1, 1)]
            / self.plasma_solver.ion_number_density.loc[(1, 1)]
        ).values
        print("Ion Ratio:", ion_ratio, ion_ratio**-1)
        print("Plasma Ion Ratio", self.plasma_solver.ion_ratio)
        ion_ratio_conv = (
            np.fabs(self.plasma_solver.ion_ratio - ion_ratio**-1)
            / ion_ratio**-1
        )
        print("Ion Ratio Conv:", ion_ratio_conv)

    def solve_continuum_state(self, continuum_estimators):
        self.base_continuum = BaseContinuum(
            plasma_array=self.plasma_solver,
            atom_data=self.atom_data,
            ws=self.simulation_state.dilution_factor,
            radiative_transition_probabilities=self.plasma_solver.transition_probabilities,
            estimators=continuum_estimators,
        )

        self.atom_data.continuum_data.set_level_number_density(
            self.plasma_solver.level_number_density
        )

        self.atom_data.continuum_data.set_level_number_density_ratio(
            self.plasma_solver
        )

    def normalize_continuum_estimators(
        self, continuum_estimators, j_blues, j_estimators
    ):
        photo_ion_norm_factor = (
            1.0
            / (
                self.transport_state.time_of_simulation
                * self.transport_state.geometry_state.volume
                * const.h.cgs.value
            ).value
        )
        damp = self.get_radiation_field_damping_factor(j_estimators)

        continuum_estimators["photo_ion_estimator"] *= (
            photo_ion_norm_factor * damp
        )
        continuum_estimators["stim_recomb_estimator"] *= (
            photo_ion_norm_factor * damp
        )
        continuum_estimators["bf_heating_estimator"] *= (
            photo_ion_norm_factor * const.h.cgs.value * damp
        )
        continuum_estimators["stim_recomb_cooling_estimator"] *= (
            photo_ion_norm_factor * const.h.cgs.value * damp
        )

        ff_norm_factor = self.get_ff_heating_norm_factor(
            self.plasma_solver.ion_number_density,
            self.plasma_solver.electron_densities.values,
            self.plasma_solver.t_electrons,
        )
        ff_norm_factor *= photo_ion_norm_factor * const.h.cgs.value * damp
        continuum_estimators["ff_heating_estimator"] *= ff_norm_factor
        j_blues *= damp
        return continuum_estimators, j_blues

    def get_radiation_field_damping_factor(self, j_estimators):
        J = (
            self.simulation_state.dilution_factor
            * self.simulation_state.t_radiative.value**4
            * const.sigma_sb.cgs.value
            / np.pi
        )
        J_estim = j_estimators / (
            4.0
            * np.pi
            * self.transport_state.time_of_simulation.value
            * self.transport_state.geometry_state.volume
        )
        damping_factor = J / J_estim
        return damping_factor

    @staticmethod
    def get_ff_heating_norm_factor(
        ion_number_density, electron_densities, t_electrons
    ):
        ionic_charge_squared = np.square(
            ion_number_density.index.get_level_values(1).values
        )
        norm_factor = (
            electron_densities
            * ion_number_density.multiply(ionic_charge_squared, axis=0)
            .sum()
            .values
        ) ** -1
        norm_factor *= np.sqrt(t_electrons)
        return norm_factor

    def solve_opacity(self):
        """Solves the opacity state and any associated objects

        Returns
        -------
        dict
            opacity_state : tardis.opacities.opacity_state.OpacityState
                State of the line opacities
            macro_atom_state : tardis.opacities.macro_atom.macro_atom_state.MacroAtomState or None
                State of the macro atom
        """
        opacity_state = self.opacity_solver.solve(self.plasma_solver)

        if self.completed_iterations == 0:
            inputs = ContinuumInputData(
                self.plasma_solver.atomic_data,
                self.plasma_solver,
                self.simulation_state.radiation_field_state.dilution_factor,
                self.plasma_solver.transition_probabilities,
                None,
            )

            photoion_rates_solver = RadiativeIonization(inputs)
            recomb_rates_solver = RadiativeRecombination(inputs)

            photoion_rate = photoion_rates_solver._calculate_rate_coefficient()
            recombination_rate = (
                recomb_rates_solver._calculate_rate_coefficient()
            )

            macro_atom_state = self.macro_atom_solver.solve(
                self.plasma_solver.j_blues,
                opacity_state.beta_sobolev,
                self.plasma_solver.stimulated_emission_factor,
                photoion_rate,
                recombination_rate,
                self.plasma_solver.coll_deexc_coeff,
                self.plasma_solver.coll_exc_coeff,
                self.plasma_solver.electron_densities,
                self.plasma_solver.level_number_density,
                self.plasma_solver.delta_E_yg,
            )

        else:
            montecarlo_globals.CONTINUUM_PROCESSES_ENABLED = True
            j_blues_df = pd.DataFrame(
                self.transport_state.estimators_line.mean_intensity_blue,
                index=self.plasma_solver.lines.index,
            )
            macro_atom_state = self.macro_atom_solver.solve(
                j_blues_df,
                opacity_state.beta_sobolev,
                self.plasma_solver.stimulated_emission_factor,
                self.plasma_solver.gamma_corr,
                self.plasma_solver.alpha_sp,
                self.plasma_solver.coll_deexc_coeff,
                self.plasma_solver.coll_exc_coeff,
                self.plasma_solver.electron_densities,
                self.plasma_solver.level_number_density,
                self.plasma_solver.delta_E_yg,
            )

        return {
            "opacity_state": opacity_state,
            "macro_atom_state": macro_atom_state,
        }

    def solve_montecarlo(
        self, opacity_states, no_of_real_packets, no_of_virtual_packets=0
    ):
        """Solve the MonteCarlo process

        Parameters
        ----------
        opacity_states : dict
            Opacity and (optionally) Macro Atom states.
        no_of_real_packets : int
            Number of real packets to simulate
        no_of_virtual_packets : int, optional
            Number of virtual packets to simulate per interaction, by default 0

        Returns
        -------
        MonteCarloTransportState
            The new transport state after simulation
        ndarray
            Array of unnormalized virtual packet energies in each frequency bin
        """
        opacity_state = opacity_states["opacity_state"]
        macro_atom_state = opacity_states["macro_atom_state"]

        self.transport_state = self.transport_solver.initialize_transport_state(
            self.simulation_state,
            opacity_state,
            macro_atom_state,
            self.plasma_solver,
            no_of_real_packets,
            no_of_virtual_packets=no_of_virtual_packets,
            iteration=self.completed_iterations,
        )

        virtual_packet_energies = self.transport_solver.run(
            self.transport_state,
            show_progress_bars=self.show_progress_bars,
        )

        output_energy = self.transport_state.packet_collection.output_energies
        if np.sum(output_energy < 0) == len(output_energy):
            logger.critical("No r-packet escaped through the outer boundary.")

        return virtual_packet_energies

    def initialize_spectrum_solver(
        self,
        opacity_states,
        virtual_packet_energies=None,
    ):
        """Set up the spectrum solver

        Parameters
        ----------
        virtual_packet_energies : ndarray, optional
            Array of virtual packet energies binned by frequency, by default None
        """
        # Set up spectrum solver
        self.spectrum_solver.transport_state = self.transport_state

        if virtual_packet_energies is not None:
            self.spectrum_solver._montecarlo_virtual_luminosity.value[:] = (
                virtual_packet_energies
            )

        if self.integrated_spectrum_settings is not None:
            # Set up spectrum solver integrator
            self.spectrum_solver.integrator_settings = (
                self.integrated_spectrum_settings
            )
            integrator_settings = self.spectrum_solver.integrator_settings
            formal_integrator = FormalIntegralSolver(
                integrator_settings.points,
                integrator_settings.interpolate_shells,
                getattr(integrator_settings, "method", None),
            )
            self.spectrum_solver.setup_optional_spectra(
                self.transport_state,
                virtual_packet_luminosity=None,
                integrator=formal_integrator,
                simulation_state=self.simulation_state,
                transport=self.transport_solver,
                plasma=self.plasma_solver,
                opacity_state=opacity_states["opacity_state"],
                macro_atom_state=opacity_states["macro_atom_state"],
            )

    def run(self):
        """Run the TARDIS simulation until convergence is reached"""
        # Initialize iterations progress bar if showing progress bars
        if self.show_progress_bars:
            initialize_iterations_pbar(self.total_iterations)

        self.converged = False

        while self.completed_iterations < self.total_iterations - 1:
            logger.info(
                f"\n\tStarting iteration {(self.completed_iterations + 1):d} of {self.total_iterations:d}"
            )

            print("Solving opacity")

            self.opacity_states = self.solve_opacity()

            print("Solving Monte Carlo transport")
            virtual_packet_energies = self.solve_montecarlo(
                self.opacity_states, self.real_packet_count
            )

            (
                estimated_values,
                estimated_radfield_properties,
            ) = self.get_convergence_estimates()

            print("Updating simulation")
            self.solve_simulation_state(estimated_values)

            continuum_estimators, j_blues = self.update_continuum_estimators()

            print("Solving plasma")
            self.solve_plasma(continuum_estimators, j_blues)

            # After first MC step
            print("Solving thermal balance")
            self.solve_thermal_balance()

            self.solve_continuum_state(continuum_estimators)

            self.converged = self.check_convergence(estimated_values)
            self.completed_iterations += 1

            if self.converged and self.convergence_strategy.stop_if_converged:
                break

        if self.converged:
            logger.info("\n\tStarting final iteration")
        else:
            logger.error(
                "\n\tITERATIONS HAVE NOT CONVERGED, starting final iteration"
            )
        self.opacity_states = self.solve_opacity()
        virtual_packet_energies = self.solve_montecarlo(
            self.opacity_states,
            self.final_iteration_packet_count,
            self.virtual_packet_count,
        )

        self.initialize_spectrum_solver(
            self.opacity_states,
            virtual_packet_energies,
        )
