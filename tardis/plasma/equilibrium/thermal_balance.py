class ThermalBalanceSolver:
    def __init__(self, thermal_rates: dict):
        self.thermal_rates = thermal_rates

    def solve(
        self,
        thermal_electron_distribution,
        radiation_field,
        level_population,
        ion_population,
        collisional_ionization_rate_coefficient,
        collisional_deexcitation_rate_coefficient,
        collisional_excitation_rate_coefficient,
        free_free_heating_estimator,
        saha_factor,
    ):
        """Compute the current heating rate and the fractional heating
        rate using all available processes. See section 6.4 in Lucy 03.

        Parameters
        ----------
        thermal_electron_distribution : ThermalElectronEnergyDistribution
            Electron energy, temperature, and density.
        radiation_field : RadiationField
            Radiation field that can calculate its own mean intensity.
        level_population : pd.DataFrame
            Level number density.
        ion_population : pd.DataFrame
            Ion number density.
        collisional_ionization_rate_coefficient : pd.DataFrame
            Collisional ionization rate coefficient.
        collisional_deexcitation_rate_coefficient : pd.DataFrame
            Collisional deexcitation rate coefficient.
        collisional_excitation_rate_coefficient : pd.DataFrame
            Collisional excitation rate coefficient.
        free_free_heating_estimator : pd.DataFrame
            Montecarlo estimator for free-free heating.
        saha_factor : pd.DataFrame
            Saha factor. See Lucy 03, equation 14.

        Returns
        -------
        pd.Series, pd.Series
            Total heating rate and fractional heating rate for each cell.
        """
        electron_density = thermal_electron_distribution.density

        bound_free_heating, free_bound_cooling = self.thermal_rates[
            "bound_free"
        ].solve(
            level_population,
            ion_population,
            thermal_electron_distribution,
            radiation_field,
            saha_factor,
        )

        free_free_heating, free_free_cooling = self.thermal_rates[
            "free_free"
        ].solve(
            free_free_heating_estimator,
            thermal_electron_distribution,
            ion_population,
        )

        collisional_ionization_heating, collisional_ionization_cooling = (
            self.thermal_rates["collisional_ionization"].solve(
                electron_density,
                ion_population,
                level_population,
                collisional_ionization_rate_coefficient,
                saha_factor,
            )
        )

        collisional_bound_heating, collisional_bound_cooling = (
            self.thermal_rates["collisional_bound"].solve(
                electron_density,
                collisional_deexcitation_rate_coefficient,
                collisional_excitation_rate_coefficient,
                level_population,
            )
        )

        total_heating = (
            bound_free_heating
            + free_free_heating
            + collisional_ionization_heating
            + collisional_bound_heating
        )
        total_cooling = (
            free_bound_cooling
            + free_free_cooling
            + collisional_ionization_cooling
            + collisional_bound_cooling
        )

        total_heating_rate = total_heating - total_cooling
        fractional_heating_rate = (
            total_heating - total_cooling
        ) / total_cooling

        return total_heating_rate, fractional_heating_rate
