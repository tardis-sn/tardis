import logging
from collections.abc import Mapping

import numpy as np
import numpy.typing as npt
import pandas as pd
from astropy import units as u
from scipy.interpolate import interp1d
from scipy.optimize import least_squares as lsq
from scipy.sparse import block_diag

from tardis import constants as const
from tardis.io.atom_data.parse_atom_data import parse_atom_data
from tardis.model import SimulationState
from tardis.opacities.continuum.continuum_state import ContinuumOpacityState
from tardis.opacities.continuum.macro_atom_state import (
    ContinuumMacroAtomState,
)
from tardis.opacities.macro_atom.macroatom_solver import (
    ContinuumMacroAtomSolver,
)
from tardis.opacities.opacity_solver import OpacitySolver
from tardis.opacities.tau_sobolev import (
    SOBOLEV_COEFFICIENT,
    calculate_beta_sobolev,
    calculate_sobolev_line_opacity,
)
from tardis.plasma.assembly import PlasmaSolverFactory
from tardis.plasma.equilibrium.evaluator import (
    PlasmaEquilibriumEvaluation,
    PlasmaEquilibriumEvaluator,
    calculate_lte_populations,
)
from tardis.plasma.equilibrium.inputs import (
    ContinuumCoefficientState,
    NumberDensityPerShell,
    SobolevInputs,
)
from tardis.plasma.equilibrium.ion_populations import IonPopulationSolver
from tardis.plasma.equilibrium.rate_matrix import (
    AnalyticIonRateMatrix,
    EstimatedIonRateMatrix,
    RateMatrix,
)
from tardis.plasma.equilibrium.rates import (
    AnalyticPhotoionizationRateSolver,
    CollisionalIonizationRateSolver,
    EstimatedPhotoionizationRateSolver,
    ThermalCollisionalRateSolver,
)
from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    BoundFreeThermalRates,
    CollisionalBoundThermalRates,
    CollisionalIonizationThermalRates,
    FreeFreeThermalRates,
)
from tardis.plasma.equilibrium.rates.radiative_rates import RadiativeRatesSolver
from tardis.plasma.equilibrium.thermal_balance import ThermalBalanceSolver
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.simulation.convergence import ConvergenceSolver
from tardis.spectrum.base import SpectrumSolver
from tardis.spectrum.luminosity import (
    calculate_filtered_luminosity,
)
from tardis.transport.montecarlo.configuration import montecarlo_globals
from tardis.transport.montecarlo.estimators import init_estimators_continuum
from tardis.transport.montecarlo.modes.iip.solver import (
    MCTransportSolverIIP,
)
from tardis.transport.montecarlo.progress_bars import initialize_iterations_pbar
from tardis.util.environment import Environment
from tardis.workflows.workflow_logger import WorkflowLogger

# logging support
logger = logging.getLogger(__name__)


class TypeIIPWorkflow:
    """Coordinate Type IIP plasma, transport, and thermal-balance solves."""

    show_progress_bars = Environment.allows_widget_display()
    enable_virtual_packet_logging = False
    log_level = None
    specific_log_level = None

    def __init__(self, configuration, csvy=False):
        """Initialize a TARDIS workflow for Type IIP supernovae.

        Parameters
        ----------
        configuration : Configuration
            Configuration object for the simulation
        csvy : bool, optional
            Set true if the configuration uses CSVY, by default False
        """
        self.workflow_logger = WorkflowLogger(
            configuration, self.log_level, self.specific_log_level
        )
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

        self.configuration = configuration

        montecarlo_globals.CONTINUUM_PROCESSES_ENABLED = True

        continuum_interactions = configuration.plasma.continuum_interaction
        configured_continuum_species = {
            tuple(species) if isinstance(species, tuple) else species
            for species in continuum_interactions.species
        }
        if not {"H I", (1, 0)} & configured_continuum_species:
            raise ValueError(
                "Continuum interactions for 'H I' must be included for the "
                "IIP workflow. Check plasma.continuum_interaction.species in "
                "the configuration."
            )

        elemental_number_density = (
            self.simulation_state.calculate_elemental_number_density(
                self.atom_data.atom_data.mass
            )
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

        factory = PlasmaSolverFactory(self.atom_data, configuration)
        factory.helium_treatment = "none"
        factory.legacy_nlte_species = list(
            dict.fromkeys([*factory.legacy_nlte_species, (1, 0)])
        )
        factory.prepare_factory(
            self.simulation_state.abundance.index,
            "tardis.plasma.properties.property_collections",
            configuration,
            allow_continuum=True,
        )
        if (1, 0, 0) in self.atom_data.photoionization_data.index:
            self.atom_data.photoionization_data.loc[(1, 0, 0), "x_sect"] = 0.0
        self.plasma_solver = factory.assemble(
            elemental_number_density,
            radiation_field,
            configuration.supernova.time_explosion.to("s"),
            link_t_rad_t_electron=(
                configuration.plasma.link_t_rad_t_electron
                * np.ones(self.simulation_state.geometry.no_of_shells_active)
            ),
        )
        self.plasma_solver.freeze(
            "electron_densities",
            "ion_number_density",
            "level_number_density",
        )
        self._continuum_estimators = None
        self._tau_sobolev = calculate_sobolev_line_opacity(
            self.atom_data.lines,
            self.plasma_solver.level_number_density,
            self.plasma_solver.time_explosion,
            self.plasma_solver.stimulated_emission_factor,
        )
        self._beta_sobolev = calculate_beta_sobolev(self._tau_sobolev)
        maximum_electron_density = (
            self.plasma_solver.number_density.multiply(
                self.plasma_solver.number_density.index.values, axis=0
            )
            .sum()
            .to_numpy()
        )
        initial_evaluator = self._build_thermal_balance_evaluator(
            maximum_electron_density, analytic=True
        )
        calculated_continuum_coefficients = (
            initial_evaluator.calculate_continuum_coefficients(
                self.plasma_solver.t_electrons
            )
        )
        initial_continuum_coefficients = calculated_continuum_coefficients[4]
        initial_level_to_continuum_saha_factor = (
            calculated_continuum_coefficients[1]
        )
        self._build_continuum_states(
            initial_continuum_coefficients,
            initial_level_to_continuum_saha_factor,
        )

        line_interaction_type = configuration.plasma.line_interaction_type

        self.opacity_solver = OpacitySolver(
            line_interaction_type,
            configuration.plasma.disable_line_scattering,
        )

        self.macro_atom_solver = ContinuumMacroAtomSolver(
            self.atom_data.levels,
            self.atom_data.lines,
            self.atom_data.photoionization_data,
            self.atom_data.ionization_data,
            selected_continuum_transitions=np.asarray(
                factory.continuum_interaction_species_multi_index.tolist()
            ),
            line_interaction_type=line_interaction_type,
            nthreads=configuration.montecarlo.nthreads,
        )

        self.transport_state = None
        self.transport_solver = MCTransportSolverIIP.from_config(
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

    def _build_continuum_states(
        self,
        continuum_coefficients: ContinuumCoefficientState,
        level_to_continuum_saha_factor: pd.DataFrame,
    ) -> None:
        """Build macro-atom and transport state from accepted equilibrium."""
        plasma = self.plasma_solver
        self._continuum_coefficients = continuum_coefficients
        self._level_to_continuum_saha_factor = level_to_continuum_saha_factor
        stimulated_cooling = None
        if self._continuum_estimators is not None:
            stimulated_cooling = pd.DataFrame(
                np.asarray(
                    self._continuum_estimators["stim_recomb_cooling_estimator"]
                ),
                index=plasma.level2continuum_idx.index,
                columns=plasma.level_number_density.columns,
            )
        self.continuum_macro_atom_state = (
            ContinuumMacroAtomState.from_equilibrium(
                plasma.atomic_data,
                plasma.lines,
                plasma.photo_ion_cross_sections,
                continuum_coefficients,
                level_to_continuum_saha_factor,
                plasma.ion_number_density,
                plasma.level_number_density,
                plasma.electron_densities,
                plasma.t_electrons,
                stimulated_cooling,
            )
        )
        self.continuum_opacity_state = ContinuumOpacityState.from_equilibrium(
            plasma.nu_i,
            plasma.level2continuum_idx,
            plasma.photo_ion_cross_sections,
            plasma.photo_ion_idx,
            plasma.ion_number_density,
            plasma.level_number_density,
            plasma.electron_densities,
            plasma.t_electrons,
            level_to_continuum_saha_factor,
            self.continuum_macro_atom_state.radiative_recombination_rate,
        )

    @staticmethod
    def initialize_radiation_field(
        geometry_state, number_density, initial_t_inner, dilution_factor
    ):
        """Set up the radiation field properties for a IIP

        Parameters
        ----------
        geometry_state : HomologousRadial1DGeometry
            The geometry state of the simulation
        number_density : pd.DataFrame
            Elemental number density
        initial_t_inner : float
            Initial inner boundary temperature
        dilution_factor : np.ndarray
            Initial dilution factor

        Returns
        -------
        np.ndarray, np.ndarray
            Radiative temperature, dilution factor arrays for the radiation field
        """
        r_inner = geometry_state.r_inner_active.value
        r_outer = geometry_state.r_outer_active.value
        r_middle = (
            geometry_state.r_middle_active.value
        )  # (r_outer + r_inner) / 2.0
        delta_r = r_outer - r_inner
        t_inner = initial_t_inner

        sigma_T = const.sigma_T.cgs.value

        N_H = number_density.loc[1].values

        # alternative tau calculation from ctardis
        # v_inner = geometry_state.v_inner_active.value
        # doppler_factor = 1.0 - v_inner / const.c.cgs.value
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
        dilution_factor[flat_mask] = 1.0
        lin_mask = np.logical_and(
            np.logical_not(flat_mask), np.logical_not(geometric_mask)
        )
        dilution_factor[lin_mask] = a * np.log(tau_middle[lin_mask]) + b

        return t_rads, dilution_factor

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
                self.transport_state.geometry_state_numba.volume,
                self.transport_state.opacity_state_numba.line_list_nu,
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

        logger.info(
            "\n\tLuminosity emitted   = %.3e\n\tLuminosity absorbed  = %.3e\n",
            emitted_luminosity,
            absorbed_luminosity,
        )

        self.workflow_logger.log_plasma_state(
            self.simulation_state.t_radiative,
            self.simulation_state.dilution_factor,
            self.simulation_state.t_inner,
            estimated_t_radiative,
            estimated_dilution_factor,
            self.simulation_state.t_inner,
        )
        # ctardis does not update t_inner
        return {
            "t_radiative": estimated_t_radiative,
            "dilution_factor": estimated_dilution_factor,
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
            no_of_shells = self.simulation_state.no_of_shells

            convergence_statuses.append(
                solver.get_convergence_status(
                    current_value, estimated_value, no_of_shells
                )
            )

        if np.all(convergence_statuses):
            hold_iterations = self.convergence_strategy.hold_iterations
            self.consecutive_converges_count += 1
            logger.info(
                "Iteration converged %d/%d consecutive times.",
                self.consecutive_converges_count,
                hold_iterations + 1,
            )
            return self.consecutive_converges_count >= hold_iterations + 1

        self.consecutive_converges_count = 0
        return False

    def solve_simulation_state(
        self,
        estimated_values,
    ):
        """Update the simulation state from the previous iteration.

        Apply the configured convergence solvers to the new estimates.

        Parameters
        ----------
        estimated_values : dict
            Estimated from the previous iterations

        Returns
        -------
        dict
            Updated values for the simulation state
        """
        next_values = {}

        for estimate_name, solver in self.convergence_solvers.items():
            next_values[estimate_name] = solver.converge(
                getattr(self.simulation_state, estimate_name),
                estimated_values[estimate_name],
            )

        self.simulation_state.t_radiative = next_values["t_radiative"]
        self.simulation_state.dilution_factor = next_values["dilution_factor"]

        return next_values

    def update_estimators(self):
        """Update the estimators for the radiation field

        Returns
        -------
        dict, pd.DataFrame
            Continuum interaction estimators and J_blues DataFrame for the radiation field
        """
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

        damped_j_blues = (
            self.transport_solver.radfield_prop_solver.estimate_jblues(
                self.transport_state.estimators_line.mean_intensity_blueward,
                self.simulation_state.radiation_field_state,
                self.transport_state.time_explosion,
                self.transport_state.time_of_simulation,
                self.transport_state.geometry_state_numba.volume,
                self.transport_state.opacity_state_numba.line_list_nu,
                detailed_optical_window=True,
            )
        )

        damped_j_blues_df = pd.DataFrame(
            damped_j_blues,
            index=self.plasma_solver.lines.index,
        )

        continuum_estimators, damped_normalized_j_blues_df = (
            self.normalize_continuum_estimators(
                continuum_estimators,
                damped_j_blues_df,
                self.transport_state.estimators_bulk.mean_intensity_total,
            )
        )

        return continuum_estimators, damped_normalized_j_blues_df

    def solve_plasma(
        self,
        continuum_estimators: Mapping[str, object],
        j_blues_df: pd.DataFrame,
    ) -> None:
        """Update the plasma solution with the new radiation field estimates

        Parameters
        ----------
        continuum_estimators : dict
            Continuum interaction estimators
        j_blues_df : pd.DataFrame
            J_blues DataFrame for the radiation field
        """
        radiation_field = DilutePlanckianRadiationField(
            self.simulation_state.t_radiative,
            self.simulation_state.dilution_factor,
        )
        self.simulation_state.radiation_field_state = radiation_field
        self._continuum_estimators = continuum_estimators
        self.plasma_solver.update(
            dilute_planckian_radiation_field=radiation_field,
            j_blues=j_blues_df,
        )

    def _build_thermal_balance_evaluator(
        self,
        maximum_electron_density: npt.NDArray[np.float64],
        *,
        analytic: bool = False,
        fixed_estimators: object | None = None,
    ) -> PlasmaEquilibriumEvaluator:
        """Build the fixed evaluator snapshot for one outer thermal solve."""
        plasma = self.plasma_solver

        if (
            not analytic
            and fixed_estimators is None
            and self._continuum_estimators is None
        ):
            raise RuntimeError(
                "Thermal balance requires a completed Monte Carlo estimator "
                "snapshot."
            )

        # The workflow has already normalized the estimator arrays. Unit time
        # and volume preserve those coefficients while adapting them to the
        # standard estimator container expected by the evaluator.
        time_simulation = 1.0 * u.s
        volume = 1.0 * u.cm**3
        estimator_scale = const.h.cgs.value
        continuum_estimators = self._continuum_estimators
        estimators = fixed_estimators
        if not analytic and estimators is None:
            estimators = init_estimators_continuum(
                np.asarray(continuum_estimators["photo_ion_estimator"]).shape,
                len(plasma.number_density.columns),
            )
            estimators.photo_ion_estimator[:] = (
                np.asarray(continuum_estimators["photo_ion_estimator"])
                * estimator_scale
            )
            estimators.stim_recomb_estimator[:] = (
                np.asarray(continuum_estimators["stim_recomb_estimator"])
                * estimator_scale
            )
            estimators.bf_heating_estimator[:] = np.asarray(
                continuum_estimators["bf_heating_estimator"]
            )
            estimators.stim_recomb_cooling_estimator[:] = np.asarray(
                continuum_estimators["stim_recomb_cooling_estimator"]
            )
            estimators.ff_heating_estimator[:] = np.asarray(
                continuum_estimators["ff_heating_estimator"]
            )

        level2continuum_edge_idx = plasma.level2continuum_idx
        photoionization_data = plasma.photo_ion_cross_sections
        equilibrium_levels = plasma.atomic_data.levels.loc[
            plasma.level_number_density.index
        ]
        level_index = plasma.level_number_density.index
        hydrogen_species = (1, 0)
        hydrogen_level_positions = np.flatnonzero(
            (
                level_index.get_level_values("atomic_number")
                == hydrogen_species[0]
            )
            & (
                level_index.get_level_values("ion_number")
                == hydrogen_species[1]
            )
        )
        population_geometries = tuple(
            NumberDensityPerShell(
                plasma.number_density.loc[1, shell],
                plasma.level_number_density[shell].to_numpy(dtype=np.float64),
                hydrogen_level_positions,
            )
            for shell in plasma.number_density.columns
        )

        line_index = plasma.lines.index
        line_species_index = line_index.droplevel(
            ["level_number_lower", "level_number_upper"]
        )
        nlte_lines_mask = np.asarray(
            line_species_index.isin([hydrogen_species]), dtype=bool
        )
        time_explosion_seconds = plasma.time_explosion
        if isinstance(time_explosion_seconds, u.Quantity):
            time_explosion_seconds = time_explosion_seconds.to_value("s")
        sobolev_input = SobolevInputs(
            plasma.lines_lower_level_index,
            plasma.lines_upper_level_index,
            plasma.g.iloc[plasma.lines_lower_level_index].to_numpy(),
            plasma.g.iloc[plasma.lines_upper_level_index].to_numpy(),
            plasma.metastability.iloc[
                plasma.lines_upper_level_index
            ].to_numpy(),
            nlte_lines_mask,
            (
                plasma.lines.wavelength_cm.to_numpy()
                * plasma.lines.f_lu.to_numpy()
                * SOBOLEV_COEFFICIENT
                * time_explosion_seconds
            ),
            np.arange(len(line_index), dtype=np.int64),
            line_index,
        )
        rate_matrix_solver = RateMatrix(
            RadiativeRatesSolver(plasma.lines),
            ThermalCollisionalRateSolver(
                equilibrium_levels,
                plasma.lines,
                plasma.atomic_data.collision_data_temperatures,
                plasma.atomic_data.yg_data,
                collision_strengths_type="cmfgen",
            ),
            equilibrium_levels,
        )
        lte_ion_population, lte_level_population = calculate_lte_populations(
            plasma.thermal_phi_lte,
            plasma.thermal_lte_partition_function,
            plasma.number_density,
            plasma.electron_densities,
            plasma.thermal_lte_level_boltzmann_factor,
            equilibrium_levels,
        )
        if analytic:
            ion_rate_matrix = AnalyticIonRateMatrix(
                AnalyticPhotoionizationRateSolver(photoionization_data),
                CollisionalIonizationRateSolver(photoionization_data),
            )
        else:
            ion_rate_matrix = EstimatedIonRateMatrix(
                EstimatedPhotoionizationRateSolver(
                    photoionization_data,
                    level2continuum_edge_idx,
                    estimators,
                    time_simulation,
                    volume,
                ),
                CollisionalIonizationRateSolver(photoionization_data),
                plasma.phi,
            )
        ion_population_solver = IonPopulationSolver(ion_rate_matrix)
        collisional_index = rate_matrix_solver.electron_rate_solver.all_collisional_strengths_index
        lower_level_index = pd.MultiIndex.from_arrays(
            [
                collisional_index.get_level_values("atomic_number"),
                collisional_index.get_level_values("ion_number"),
                collisional_index.get_level_values("level_number_lower"),
            ],
            names=["atomic_number", "ion_number", "level_number"],
        )
        upper_level_index = pd.MultiIndex.from_arrays(
            [
                collisional_index.get_level_values("atomic_number"),
                collisional_index.get_level_values("ion_number"),
                collisional_index.get_level_values("level_number_upper"),
            ],
            names=["atomic_number", "ion_number", "level_number"],
        )
        collisional_line_frequencies = (
            equilibrium_levels.energy.loc[upper_level_index].to_numpy()
            - equilibrium_levels.energy.loc[lower_level_index].to_numpy()
        ) / const.h.cgs.value
        thermal_balance_solver = None
        thermal_balance_arguments = None
        if not analytic and fixed_estimators is None:
            thermal_balance_solver = ThermalBalanceSolver(
                BoundFreeThermalRates(photoionization_data),
                FreeFreeThermalRates(),
                CollisionalIonizationThermalRates(photoionization_data),
                CollisionalBoundThermalRates(
                    pd.DataFrame({"nu": collisional_line_frequencies})
                ),
            )
            thermal_balance_arguments = {
                "free_free_heating_estimator": np.asarray(
                    continuum_estimators["ff_heating_estimator"]
                ),
                "bound_free_heating_estimator": pd.DataFrame(
                    np.asarray(continuum_estimators["bf_heating_estimator"]),
                    index=level2continuum_edge_idx.index,
                    columns=plasma.number_density.columns,
                ),
                "stimulated_recombination_estimator": pd.DataFrame(
                    np.asarray(
                        continuum_estimators["stim_recomb_cooling_estimator"]
                    ),
                    index=level2continuum_edge_idx.index,
                    columns=plasma.number_density.columns,
                ),
            }
        return PlasmaEquilibriumEvaluator(
            photoionization_data,
            level2continuum_edge_idx,
            estimators,
            time_simulation,
            volume,
            equilibrium_levels,
            plasma.ionization_data,
            rate_matrix_solver,
            pd.DataFrame(
                plasma.j_blues,
                index=line_index,
                columns=plasma.number_density.columns,
            ),
            population_geometries,
            tuple(sobolev_input for _ in plasma.number_density.columns),
            level_index,
            hydrogen_species,
            plasma.number_density,
            maximum_electron_density,
            ion_population_solver=ion_population_solver,
            ion_population_arguments={
                "radiation_field": (
                    plasma.dilute_planckian_radiation_field
                    if analytic
                    else None
                ),
                "elemental_number_density": plasma.number_density,
                "lte_level_population": lte_level_population,
                "lte_ion_population": lte_ion_population,
                "estimated_ion_population": plasma.ion_number_density,
                "partition_function": plasma.thermal_lte_partition_function,
                "boltzmann_factor": plasma.thermal_lte_level_boltzmann_factor,
            },
            thermal_balance_solver=thermal_balance_solver,
            thermal_balance_arguments=thermal_balance_arguments,
            radiation_field=(
                plasma.dilute_planckian_radiation_field if analytic else None
            ),
        )

    def _initialize_thermal_balance_evaluator(
        self,
        maximum_electron_density: npt.NDArray[np.float64],
    ) -> None:
        """Freeze the evaluator and level seed for one outer solve."""
        plasma = self.plasma_solver
        hydrogen_species = (1, 0)
        self._thermal_balance_evaluator = self._build_thermal_balance_evaluator(
            maximum_electron_density
        )
        self._thermal_balance_radiation_temperature = np.asarray(
            plasma.t_rad, dtype=np.float64
        ).copy()
        self._thermal_balance_level_seed = plasma.level_number_density.loc[
            hydrogen_species
        ].divide(plasma.ion_number_density.loc[hydrogen_species], axis=1)

    def thermal_balance_iteration(
        self,
        initial_guess: npt.NDArray[np.float64],
        max_electron_number_density: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """Compute the fixed-candidate thermal-balance residual.

        Parameters
        ----------
        initial_guess : np.ndarray
            Initial guess for the electron number density fraction and
            link_t_rad_t_electron for each shell, in the form
            [n_e_frac_1, link_1, n_e_frac_2, link_2,...]
        max_electron_number_density : np.ndarray
            Maximum possible electron number density for each shell.

        Returns
        -------
        np.ndarray
            Interleaved normalized electron and fractional-heating residuals.
        """
        electron_density_fraction = initial_guess[::2]
        link_t_rad_t_electron = initial_guess[1::2]

        logger.info("Link: %s", link_t_rad_t_electron)

        electron_densities = (
            max_electron_number_density * electron_density_fraction
        )
        evaluation = self._thermal_balance_evaluator.evaluate(
            electron_densities,
            np.asarray(
                self._thermal_balance_radiation_temperature
                * link_t_rad_t_electron,
                dtype=np.float64,
            ),
            self._thermal_balance_level_seed,
        )

        electron_residual = evaluation.electron_residual.to_numpy()
        fractional_heating = evaluation.fractional_heating.to_numpy()
        solution = np.zeros_like(initial_guess)
        if not np.isfinite(fractional_heating).all():
            logger.warning("Heating not finite\n")
        if not np.isfinite(electron_residual).all():
            logger.warning("Fractional electron change not finite\n")

        solution[::2] = electron_residual
        solution[1::2] = fractional_heating
        logger.info(
            "Normalized electron fraction change: %s",
            electron_residual,
        )
        logger.info("Heating: %s", fractional_heating)
        return solution

    def _publish_thermal_balance_state(
        self,
        candidate: npt.NDArray[np.float64],
    ) -> None:
        """Publish one accepted evaluator result to the plasma graph."""
        evaluation = self._thermal_balance_evaluation
        self.plasma_solver.update(
            electron_densities=evaluation.charge_solved_electron_density,
            ion_number_density=evaluation.ion_population,
            level_number_density=evaluation.absolute_level_population,
            link_t_rad_t_electron=candidate[1::2],
        )
        self._tau_sobolev = evaluation.tau_sobolev
        self._beta_sobolev = evaluation.beta_sobolev

    @staticmethod
    def _validate_thermal_balance_evaluation(
        evaluation: PlasmaEquilibriumEvaluation,
    ) -> None:
        """Reject a final evaluator state that does not physically close."""
        failures = []
        population_fields = {
            "normalized_population": evaluation.normalized_population,
            "absolute_level_population": evaluation.absolute_level_population,
            "ion_population": evaluation.ion_population,
            "electron_density": evaluation.charge_solved_electron_density,
        }
        for field_name, population in population_fields.items():
            if population is None:
                failures.append(f"{field_name} is missing")
                continue
            population_values = np.asarray(population, dtype=np.float64)
            if not np.isfinite(population_values).all():
                failures.append(f"{field_name} is nonfinite")
            if (population_values < 0.0).any():
                failures.append(f"{field_name} is negative")

        normalized_totals = evaluation.normalized_population.sum(axis=0)
        if not np.allclose(normalized_totals, 1.0, rtol=1e-12, atol=0.0):
            failures.append("normalized_population does not sum to one")

        for field_name, state in {
            "tau_sobolev": evaluation.tau_sobolev,
            "beta_sobolev": evaluation.beta_sobolev,
        }.items():
            if not np.isfinite(np.asarray(state, dtype=np.float64)).all():
                failures.append(f"{field_name} is nonfinite")

        closure_tolerances = {
            "trial_level_residual": 1e-10,
            "level_residual": 1e-10,
            "charge_residual": 1e-10,
            # The outer trial density and the independently projected charge
            # solution differ at the established Phase 2 closure floor.
            "electron_residual": 2e-8,
            # Heating closes by cancellation; these measured Phase 2 floors
            # are stricter than the observable legacy-parity contract.
            "total_heating": 5e-13,
            "fractional_heating": 2e-7,
        }
        for field_name, tolerance in closure_tolerances.items():
            residual = getattr(evaluation, field_name)
            if residual is None:
                failures.append(f"{field_name} is missing")
                continue
            residual_values = np.asarray(residual, dtype=np.float64)
            if not np.isfinite(residual_values).all():
                failures.append(f"{field_name} is nonfinite")
                continue
            maximum_residual = float(np.max(np.abs(residual_values)))
            if maximum_residual > tolerance:
                failures.append(
                    f"{field_name}={maximum_residual:.3e} exceeds "
                    f"{tolerance:.3e}"
                )

        if failures:
            raise RuntimeError(
                "Thermal-balance solution failed final closure: "
                + "; ".join(failures)
            )

    def solve_thermal_balance(self) -> None:
        """Solve the plasma heating and cooling balance.

        Set the electron number density and link_t_rad_t_electron
        to values that satisfy thermal balance
        """
        link_t_rad_t_electron_start = self.plasma_solver.link_t_rad_t_electron
        if np.array_equal(
            link_t_rad_t_electron_start,
            np.ones_like(link_t_rad_t_electron_start),
        ):
            link_t_rad_t_electron_start = (
                self.simulation_state.radiation_field_state.dilution_factor
                ** 0.25
            )
            logger.info(
                "Setting initial guess for link between T_rad and T_e from dilution factor:"
            )
            logger.info(link_t_rad_t_electron_start)

        max_electron_number_density = (
            (
                self.plasma_solver.number_density.multiply(
                    self.plasma_solver.number_density.index.values, axis=0
                )
            )
            .sum()
            .values
        )
        initial_electron_fraction = (
            self.plasma_solver.electron_densities / max_electron_number_density
        ).values

        logger.info("Initial electron fraction: %s", initial_electron_fraction)

        initial_guess = np.zeros(2 * len(link_t_rad_t_electron_start))
        initial_guess[::2] = initial_electron_fraction
        initial_guess[1::2] = link_t_rad_t_electron_start
        self._initialize_thermal_balance_evaluator(max_electron_number_density)
        no_shells = self.simulation_state.geometry.no_of_shells_active

        jac_sparsity = block_diag([np.ones((2, 2))] * no_shells)
        t_floor = 1500.0 * u.K
        minimum_t_rad_link = t_floor / self.simulation_state.t_radiative.min()
        logger.info("Minimum T_rad link: %s", minimum_t_rad_link)

        lower_bound = [0.0, minimum_t_rad_link] * no_shells
        upper_bound = [1.0, 1.5] * no_shells
        lower_bound = np.asarray(lower_bound, dtype=float)
        upper_bound = np.asarray(upper_bound, dtype=float)

        out_of_bounds = (initial_guess < lower_bound) | (
            initial_guess > upper_bound
        )
        if np.any(out_of_bounds):
            offending_indices = np.flatnonzero(out_of_bounds)
            logger.warning(
                "Out-of-bounds thermal balance initial guess values; "
                "clipping indices=%s values=%s lower=%s upper=%s",
                offending_indices,
                initial_guess[out_of_bounds],
                lower_bound[out_of_bounds],
                upper_bound[out_of_bounds],
            )
            initial_guess = np.clip(initial_guess, lower_bound, upper_bound)

        thermal_lsq_result = lsq(
            self.thermal_balance_iteration,
            initial_guess,
            bounds=(lower_bound, upper_bound),
            jac_sparsity=jac_sparsity,
            xtol=1e-14,
            ftol=1e-12,
            x_scale="jac",
            verbose=1,
            max_nfev=100,
            method="trf",
            gtol=1e-14,
            args=(max_electron_number_density,),
        )
        # Preserve the frozen seed used by the optimizer for the first rebuild,
        # then use that accepted population as the canonical final-state seed.
        accepted_candidate = thermal_lsq_result.x
        accepted_seed_evaluation = self._thermal_balance_evaluator.evaluate(
            max_electron_number_density * accepted_candidate[::2],
            np.asarray(
                self._thermal_balance_radiation_temperature
                * accepted_candidate[1::2],
                dtype=np.float64,
            ),
            self._thermal_balance_level_seed,
        )
        self._thermal_balance_evaluation = (
            self._thermal_balance_evaluator.evaluate(
                max_electron_number_density * accepted_candidate[::2],
                np.asarray(
                    self._thermal_balance_radiation_temperature
                    * accepted_candidate[1::2],
                    dtype=np.float64,
                ),
                accepted_seed_evaluation.normalized_population,
            )
        )
        self._validate_thermal_balance_evaluation(
            self._thermal_balance_evaluation
        )
        self._publish_thermal_balance_state(accepted_candidate)

    def solve_continuum_state(
        self, continuum_estimators: Mapping[str, object]
    ) -> None:
        """Build continuum interaction and opacity state from accepted plasma.

        Parameters
        ----------
        continuum_estimators : dict
            Continuum interaction estimators
        """
        self._continuum_estimators = continuum_estimators
        evaluation = self._thermal_balance_evaluation
        self._build_continuum_states(
            evaluation.continuum_coefficients,
            evaluation.level_to_continuum_saha_factor,
        )

    def normalize_continuum_estimators(
        self, continuum_estimators, j_blues, j_estimators
    ):
        """Normalize the continuum estimators and J-blues.

        Apply the Monte Carlo normalization and radiation-field damping.

        Parameters
        ----------
        continuum_estimators : dict
            Continuum interaction estimators
        j_blues : pd.DataFrame
            J_blues DataFrame for the radiation field
        j_estimators : np.ndarray
            J array for the radiation field

        Returns
        -------
        dict, pd.DataFrame
            Normalized continuum interaction estimators and J_blues DataFrame
        """
        photo_ion_norm_factor = (
            1.0
            / (
                self.transport_state.time_of_simulation
                * self.transport_state.geometry_state_numba.volume
                * const.h.cgs.value
            ).value
        )
        damping_factor = self.get_radiation_field_damping_factor(j_estimators)

        continuum_estimators["photo_ion_estimator"] *= (
            photo_ion_norm_factor * damping_factor
        )
        continuum_estimators["stim_recomb_estimator"] *= (
            photo_ion_norm_factor * damping_factor
        )
        continuum_estimators["bf_heating_estimator"] *= (
            photo_ion_norm_factor * const.h.cgs.value * damping_factor
        )
        continuum_estimators["stim_recomb_cooling_estimator"] *= (
            photo_ion_norm_factor * const.h.cgs.value * damping_factor
        )

        ff_norm_factor = self.get_ff_heating_norm_factor(
            self.plasma_solver.ion_number_density,
            self.plasma_solver.electron_densities.values,
            self.plasma_solver.t_electrons,
        )
        ff_norm_factor *= (
            photo_ion_norm_factor * const.h.cgs.value * damping_factor
        )
        continuum_estimators["ff_heating_estimator"] *= ff_norm_factor
        j_blues *= damping_factor
        return continuum_estimators, j_blues

    def get_radiation_field_damping_factor(self, j_estimators):
        """Compute the radiation field damping factor

        Parameters
        ----------
        j_estimators : np.ndarray
            J array for the radiation field

        Returns
        -------
        np.ndarray
            Damping factor
        """
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
            * self.transport_state.geometry_state_numba.volume
        )
        damping_factor = J / J_estim
        return damping_factor

    @staticmethod
    def get_ff_heating_norm_factor(
        ion_number_density, electron_densities, t_electrons
    ):
        """Compute the free-free heating normalization factor

        Parameters
        ----------
        ion_number_density : pd.DataFrame
            Ion number density DataFrame from the plasma solver
        electron_densities : np.ndarray
            Electron density array from the plasma solver
        t_electrons : np.ndarray
            Electron temperature array from the plasma solver

        Returns
        -------
        np.ndarray
            Free-free heating normalization factor
        """
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
        opacity_state = self.opacity_solver.solve(
            self.plasma_solver,
            self.continuum_opacity_state,
            tau_sobolev=self._tau_sobolev,
            beta_sobolev=self._beta_sobolev,
        )
        continuum_rates = self.continuum_macro_atom_state

        macro_atom_state = self.macro_atom_solver.solve(
            self.plasma_solver.j_blues,
            opacity_state.beta_sobolev,
            self.plasma_solver.stimulated_emission_factor,
            continuum_rates,
            self.plasma_solver.electron_densities,
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

        """
        opacity_state = opacity_states["opacity_state"]
        macro_atom_state = opacity_states["macro_atom_state"]

        self.transport_state = self.transport_solver.initialize_transport_state(
            self.simulation_state,
            opacity_state,
            macro_atom_state,
            no_of_real_packets,
            no_of_virtual_packets=no_of_virtual_packets,
            iteration=self.completed_iterations,
        )

        self.transport_solver.run(
            self.transport_state,
            show_progress_bars=self.show_progress_bars,
        )

        output_energy = self.transport_state.packet_collection.output_energies
        if np.sum(output_energy < 0) == len(output_energy):
            logger.critical("No r-packet escaped through the outer boundary.")

    def initialize_spectrum_solver(
        self,
    ):
        """Set up the spectrum solver"""
        # probably needs to expand again in the future to handle formal integral
        self.spectrum_solver.transport_state = self.transport_state

    def run(self):
        """Run the TARDIS simulation until convergence is reached"""
        # Initialize iterations progress bar if showing progress bars
        if self.show_progress_bars:
            initialize_iterations_pbar(self.total_iterations)

        self.converged = False

        while self.completed_iterations < self.total_iterations - 1:
            logger.info(
                "\n\tStarting iteration %d of %d",
                self.completed_iterations + 1,
                self.total_iterations,
            )
            logger.info("Opacity solve started.")
            self.opacity_states = self.solve_opacity()
            logger.info("Opacity solve finished.")
            self.solve_montecarlo(self.opacity_states, self.real_packet_count)
            logger.info("Montecarlo solve finished.")
            (
                estimated_values,
                _estimated_radfield_properties,
            ) = self.get_convergence_estimates()

            self.solve_simulation_state(estimated_values)

            normalized_continuum_estimators, damped_normalized_j_blues = (
                self.update_estimators()
            )
            self.damped_normalized_j_blues = damped_normalized_j_blues
            self.solve_plasma(
                normalized_continuum_estimators, damped_normalized_j_blues
            )

            # After first MC step
            self.solve_thermal_balance()
            logger.info("Thermal balance solve finished.")

            self.solve_continuum_state(normalized_continuum_estimators)

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
        self.solve_montecarlo(
            self.opacity_states,
            self.final_iteration_packet_count,
            self.virtual_packet_count,
        )

        self.initialize_spectrum_solver()
