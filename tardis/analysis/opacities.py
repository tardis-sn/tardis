"""This module provides an opacity calculator class with which the opacities
and optical depth information may be extracted from Tardis runs."""
import logging
import numpy as np
import astropy.units as units
from tardis import constants as const
from astropy.modeling.models import Blackbody

logger = logging.getLogger(__name__)


class opacity_calculator(object):
    """Basic Tardis opacity and optical depth calculator

    Given the model object of a Tardis run and a frequency grid, detailed
    opacity and optical depth information may be extracted. In particular, the
    following quantities may be calculated:
        * electron-scattering opacities in each cell
        * bound-bound scattering opacities based on the expansion opacity
          formalism of Blinnikov et al. 1998 in each cell on the provided
          frequency grid
        * total opacities (i.e. electron scattering and line opacity) per each
          cell and frequency bin
        * Planck-mean total opacities
        * Planck-mean total optical depths of the model

    All these quantities are calculated on the fly once the corresponding class
    attributes are accesses. A simple caching scheme is implemented which
    stores all opacity and optical depth quantities after they have been
    calculated for the first time. Only when the base model or the parameters
    of the frequency grid change, the opacity and optical depth are
    recalculated.

    Parameters
    ----------
    mdl : tardis.model.Radial1DModel
        model object of the Tardis run
    nbins : int
        number of bins of the frequency grid (default 300)
    lam_min : astropy.units.Quantity
        lower wavelength boundary (i.e. upper frequency boundary) of the
        frequency grid (default 100 AA)
    lam_max : astropy.units.Quantity
        upper wavelength boundary (i.e. lower frequency boundary) of the
        frequency grid (default 20000 AA)
    bin_scaling : str
        'log' for logarithmic scaling of the frequency grid, 'linear' for a
        linear one (default 'log')
    """

    def __init__(
        self,
        mdl,
        nbins=300,
        lam_min=100 * units.AA,
        lam_max=2e4 * units.AA,
        bin_scaling="log",
    ):

        # Private attributes
        self._mdl = None
        self._nbins = None
        self._lam_min = None
        self._lam_max = None
        self._bin_scaling = None

        # Private derived model attributes
        self._r_inner = None
        self._r_outer = None
        self._t_exp = None
        self._nshells = None

        # Private opacity and optical depth attributes
        self._kappa_exp = None
        self._kappa_thom = None
        self._kappa_thom_grid = None
        self._kappa_tot = None
        self._planck_kappa = None
        self._planck_delta_tau = None
        self._planck_tau = None

        # Passing the arguments
        self.mdl = mdl
        self.nbins = nbins
        self.lam_min = lam_min
        self.lam_max = lam_max
        self.bin_scaling = bin_scaling

    def _reset_opacities(self):
        """Reset all private opacity and optical depth attributes to enforce a
        re-calculation once they are accessed the next time (part of caching
        scheme). Should be called whenever a basic attribute of the calculator
        changes.
        """
        self._kappa_exp = None
        self._kappa_thom = None
        self._kappa_thom_grid = None
        self._kappa_tot = None
        self._planck_kappa = None
        self._planck_delta_tau = None
        self._planck_tau = None

    def _reset_bins(self):
        """Reset all derived private attributes associated with the frequency
        grid (part of the caching scheme). Should be called every time any
        parameters of the frequency grid are changed.
        """
        self._nu_bins = None
        self._reset_opacities()

    def _reset_model(self):
        """Reset all derived private attributes associated with the input model
        (part of the caching scheme). Should be called every time the input
        model is changed.
        """
        self._t_exp = None
        self._nshells = None
        self._r_inner = None
        self._r_outer = None

        self._reset_opacities()

    @property
    def bin_scaling(self):
        return self._bin_scaling

    @bin_scaling.setter
    def bin_scaling(self, val):
        allowed_values = ["log", "linear"]
        if val not in allowed_values:
            raise ValueError(
                "wrong bin_scaling; must be "
                "among {:s}".format(",".join(allowed_values))
            )
        self._reset_bins()
        self._bin_scaling = val

    @property
    def lam_min(self):
        return self._lam_min

    @lam_min.setter
    def lam_min(self, val):
        self._reset_bins()
        try:
            val.to("AA")
        except AttributeError:
            logger.warning("lam_min provided without units; assuming AA")
            val *= units.AA
        self._lam_min = val

    @property
    def lam_max(self):
        return self._lam_max

    @lam_max.setter
    def lam_max(self, val):
        self._reset_bins()
        try:
            val.to("AA")
        except AttributeError:
            logger.warning("lam_max provided without units; assuming AA")
            val *= units.AA
        self._lam_max = val

    @property
    def mdl(self):
        return self._mdl

    @mdl.setter
    def mdl(self, val):
        self._reset_model()
        self._mdl = val

    @property
    def nshells(self):
        """number of radial shells in the model"""
        if self._nshells is None:
            self._nshells = self.mdl.model.no_of_shells
        return self._nshells

    @property
    def t_exp(self):
        """time since explosion of the model"""
        if self._t_exp is None:
            self._t_exp = self.mdl.model.time_explosion
        return self._t_exp

    @property
    def r_inner(self):
        """inner radius of the model shells"""
        if self._r_inner is None:
            self._r_inner = self.mdl.model.r_inner
        return self._r_inner

    @property
    def r_outer(self):
        """outer radius of the model shell"""
        if self._r_outer is None:
            self._r_outer = self.mdl.model.r_outer
        return self._r_outer

    @property
    def nu_bins(self):
        """frequency grid for the opacity evaluation"""
        if self._nu_bins is None:
            nu_max = self.lam_min.to("Hz", equivalencies=units.spectral())
            nu_min = self.lam_max.to("Hz", equivalencies=units.spectral())
            if self.bin_scaling == "log":
                nu_bins = (
                    np.logspace(
                        np.log10(nu_min.value),
                        np.log10(nu_max.value),
                        self.nbins + 1,
                    )
                    * units.Hz
                )
            elif self.bin_scaling == "linear":
                nu_bins = np.linspace(nu_min, nu_max, self.nbins + 1)
            self._nu_bins = nu_bins
        return self._nu_bins

    @property
    def kappa_exp(self):
        """bound-bound opacity according to the expansion opacity formalism per
        frequency bin and cell"""
        if self._kappa_exp is None:
            self._kappa_exp = self._calc_expansion_opacity()
        return self._kappa_exp

    @property
    def kappa_thom(self):
        """Thomson scattering opacity per cell"""
        if self._kappa_thom is None:
            self._kappa_thom = self._calc_thomson_scattering_opacity()
        return self._kappa_thom

    @property
    def kappa_thom_grid(self):
        """Thomson scattering opacity per frequency bin and cell"""
        if self._kappa_thom_grid is None:
            kappa_thom_grid = np.zeros((self.nbins, self.nshells)) / units.cm
            for i in range(self.nbins):
                kappa_thom_grid[i, :] = self.kappa_thom
            self._kappa_thom_grid = kappa_thom_grid
        return self._kappa_thom_grid

    @property
    def kappa_tot(self):
        """total opacity frequency bin and cell"""
        if self._kappa_tot is None:
            kappa_tot = self.kappa_exp + self.kappa_thom_grid
            self._kappa_tot = kappa_tot
        return self._kappa_tot

    @property
    def planck_kappa(self):
        """Planck-mean opacity per cell"""
        if self._planck_kappa is None:
            planck_kappa = self._calc_planck_mean_opacity()
            self._planck_kappa = planck_kappa
        return self._planck_kappa

    @property
    def planck_delta_tau(self):
        """Planck-mean optical depth of each cell"""
        if self._planck_delta_tau is None:
            planck_delta_tau = self._calc_planck_optical_depth()
            self._planck_delta_tau = planck_delta_tau
        return self._planck_delta_tau

    @property
    def planck_tau(self):
        """Planck-mean optical depth, integrated from the surface to
        corresponding inner shell radius"""
        if self._planck_tau is None:
            planck_tau = self._calc_integrated_planck_optical_depth()
            self._planck_tau = planck_tau
        return self._planck_tau

    def _calc_expansion_opacity(self):
        """Calculate the bound-bound opacity on the provided frequency grid in
        each cell.

        The expansion opacity formalism, in the particular version of Blinnikov
        et al. 1998, is used. In the supernova ejecta case (assuming perfect
        homologous expansion), the formula for the expansion opacity in the
        interval $[\nu, \nu+\Delta \nu]$ simplifies to
        \[
          \chi_{\mathrm{exp}} = \frac{\nu}{\Delta \nu} \frac{1}{c t} \sum_j
          \left(1 - \exp(-\tau_{\mathrm{S},j})\right)
        \]
        The summation involves all lines in the frequency bin.

        Returns
        -------
        kappa_exp : astropy.units.Quantity ndarray
            expansion opacity array (shape Nbins x Nshells)
        """

        index = self.mdl.plasma.tau_sobolevs.index
        line_waves = self.mdl.plasma.atomic_data.lines.loc[index]
        line_waves = line_waves.wavelength.values * units.AA

        kappa_exp = np.zeros((self.nbins, self.nshells)) / units.cm

        for i in range(self.nbins):

            lam_low = self.nu_bins[i + 1].to(
                "AA", equivalencies=units.spectral()
            )
            lam_up = self.nu_bins[i].to("AA", equivalencies=units.spectral())

            mask = np.argwhere(
                (line_waves > lam_low) * (line_waves < lam_up)
            ).ravel()
            taus = self.mdl.plasma.tau_sobolevs.iloc[mask]
            tmp = np.sum(1 - np.exp(-taus)).values
            kappa_exp[i, :] = (
                tmp
                * self.nu_bins[i]
                / (self.nu_bins[i + 1] - self.nu_bins[i])
                / (const.c * self.t_exp)
            )

        return kappa_exp.to("1/cm")

    def _calc_thomson_scattering_opacity(self):
        """Calculate the Thomson scattering opacity for each grid cell

        \[
          \chi_{\mathrm{Thomson}} = n_{\mathrm{e}} \sigma_{\mathrm{T}}
        \]

        Returns
        -------
        kappa_thom : astropy.units.Quantity ndarray
            Thomson scattering opacity (shape Nshells)
        """

        try:
            sigma_T = const.sigma_T
        except AttributeError:
            logger.warning("using astropy < 1.1.1: setting sigma_T manually")
            sigma_T = 6.65245873e-29 * units.m**2

        edens = self.mdl.plasma.electron_densities.values

        try:
            edens.to("1/cm^3")
        except AttributeError:
            logger.info("setting electron density units by hand (cm^-3)")
            edens = edens / units.cm**3

        kappa_thom = edens * sigma_T

        return kappa_thom.to("1/cm")

    def _calc_planck_mean_opacity(self):
        """Calculate the Planck-mean total opacity for each grid cell.

        See, e.g. Mihalas and Mihalas 1984, for details on the averaging
        process.

        Returns
        -------
        kappa_planck_mean : astropy.units.Quantity ndarray
            Planck-mean opacity (shape Nshells)
        """

        kappa_planck_mean = np.zeros(self.nshells) / units.cm

        for i in range(self.nshells):
            delta_nu = self.nu_bins[1:] - self.nu_bins[:-1]
            temperature = self.mdl.plasma.t_rad[i]
            bb_nu = Blackbody(temperature)

            tmp = (
                bb_nu(self.nu_bins[:-1]) * delta_nu * self.kappa_tot[:, 0]
            ).sum()
            tmp /= (bb_nu(self.nu_bins[:-1], temperature) * delta_nu).sum()

            kappa_planck_mean[i] = tmp

        return kappa_planck_mean.to("1/cm")

    def _calc_planck_optical_depth(self):
        """Calculate the Planck-mean optical depth of each grid cell

        Returns
        -------
        delta_tau : astropy.units.Quantity ndarray (dimensionless)
            Planck-mean optical depth (shape Nshells)
        """

        delta_r = self.r_outer - self.r_inner
        delta_tau = delta_r * self.planck_kappa

        return delta_tau.to("")

    def _calc_integrated_planck_optical_depth(self):
        """Calculate integrated Planck-mean optical depth

        For each cell, the optical depth integral from the inner shell radius
        to the surface is determined.

        Returns
        -------
        tau : numpy.ndarray
            integrated Planck-mean optical depth
        """
        tau = np.zeros(self.nshells)

        tau[-1] = self.planck_delta_tau[-1]
        for i in range(self.nshells - 2, -1, -1):
            tau[i] = tau[i + 1] + self.planck_delta_tau[i]

        return tau
