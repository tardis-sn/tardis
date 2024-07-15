import warnings
import numpy as np
from astropy import units as u
from tardis.io.util import HDFWriterMixin


class TARDISSpectrum(HDFWriterMixin):
    """
    TARDISSpectrum(_frequency, luminosity)

    _frequency: astropy.units.Quantity with unit 'Hz' or a length
        These are bin edges of frequency or wavelenght bins for the spectrum.

    luminosity: astropy.units.Quantity with unit Energy per second
        The luminosity in each bin of the spectrum.

    After manually adding a distance attribute, the properties 'flux_nu' and
    'flux_lambda' become available
    """

    hdf_properties = [
        "_frequency",
        "luminosity",
        "delta_frequency",
        "wavelength",
        "luminosity_density_lambda",
    ]

    def __init__(self, _frequency, luminosity):

        # Check for correct inputs
        if not _frequency.shape[0] == luminosity.shape[0] + 1:
            raise ValueError(
                "shape of '_frequency' and 'luminosity' are not compatible"
                f": '{_frequency.shape[0]}' and '{luminosity.shape[0]}'"
            )
        self._frequency = _frequency.to("Hz", u.spectral())
        self.luminosity = luminosity.to("erg / s")

        l_nu_unit = u.def_unit("erg\ s^{-1}\ Hz^{-1}", u.Unit("erg/(s Hz)"))
        l_lambda_unit = u.def_unit(
            "erg\ s^{-1}\ \\AA^{-1}", u.Unit("erg/(s AA)")
        )

        self.frequency = self._frequency[:-1]
        self.delta_frequency = self._frequency[1] - self._frequency[0]
        self.wavelength = self.frequency.to("angstrom", u.spectral())

        self.luminosity_density_nu = (
            self.luminosity / self.delta_frequency
        ).to(l_nu_unit)
        self.luminosity_density_lambda = self.f_nu_to_f_lambda(
            self.luminosity_density_nu,
        ).to(l_lambda_unit)

    @property
    def flux_nu(self):
        warnings.simplefilter("always", DeprecationWarning)
        warnings.warn(
            "TARDISSpectrum.flux_nu is deprecated, "
            "please use TARDISSpectrum.luminosity_to_flux() in the "
            "future.",
            category=DeprecationWarning,
            stacklevel=2,
        )
        warnings.simplefilter("default", DeprecationWarning)
        try:
            return self.luminosity_to_flux(
                self.luminosity_density_nu, self.distance
            )
        except AttributeError:
            flux = "flux_nu"
            raise AttributeError(
                "distance is required as attribute of"
                f'{self.__class__.__name__} to calculate "{flux}"'
            )

    @property
    def flux_lambda(self):
        warnings.simplefilter("always", DeprecationWarning)
        warnings.warn(
            "TARDISSpectrum.flux_lambda is deprecated, "
            "please use TARDISSpectrum.luminosity_to_flux() in the "
            "future.",
            category=DeprecationWarning,
            stacklevel=2,
        )
        warnings.simplefilter("default", DeprecationWarning)
        try:
            return self.luminosity_to_flux(
                self.luminosity_density_lambda, self.distance
            )
        except AttributeError:
            flux_lambda = "flux_lambda"
            raise AttributeError(
                "distance is required as attribute of"
                f'{self.__class__.__name__} to calculate "{flux_lambda}"'
            )

    @staticmethod
    def luminosity_to_flux(luminosity, distance):
        return luminosity / (4 * np.pi * distance.to("cm") ** 2)

    def f_nu_to_f_lambda(self, f_nu):
        return f_nu * self.frequency / self.wavelength

    def plot(self, ax=None, mode="wavelength", **kwargs):
        """
        Plot the Spectrum

        Parameters
        ----------
        ax : matplotlib.axes._subplots.AxesSubplot or None, optional
            Axis on which to create plot. Default value is None which will
            create plot on a new figure's axis.
        mode : {'wavelength', 'frequency'}, optional
            Quantity to plot spectrum against, either 'wavelength' or
            'frequency'. Default value is 'wavelength'.
        **kwargs : dict, optional
            ``matplotlib.lines.Line2D`` properties to style the plot: refer to
            `matplotlib documentation
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D>`_
            for a list of all possible arguments.
        """
        if ax is None:
            from matplotlib.pyplot import gca

            ax = gca()
        if mode == "wavelength":
            ax.plot(
                self.wavelength.value,
                self.luminosity_density_lambda.value,
                **kwargs,
            )
            ax.set_xlabel(
                f"Wavelength [{self.wavelength.unit.to_string('latex_inline')}]"
            )
            ax.set_ylabel(
                f"$L_\\lambda$ [{self.luminosity_density_lambda.unit.to_string('latex_inline')}]"
            )
        elif mode == "frequency":
            ax.plot(
                self.frequency.value, self.luminosity_density_nu.value, **kwargs
            )
            ax.set_xlabel(
                f"Frequency [{self.frequency.unit.to_string('latex_inline')}]"
            )
            ax.set_ylabel(
                f"$L_\\nu$ [{self.luminosity_density_nu.unit.to_string('latex_inline')}]"
            )
        else:
            warnings.warn(f"Did not find plotting mode {mode}, doing nothing.")

    def to_ascii(self, fname, mode="luminosity_density"):
        if mode == "luminosity_density":
            np.savetxt(
                fname,
                list(
                    zip(
                        self.wavelength.value,
                        self.luminosity_density_lambda.value,
                    )
                ),
            )
        elif mode == "flux":
            np.savetxt(
                fname, zip(self.wavelength.value, self.flux_lambda.value)
            )
        else:
            raise NotImplementedError(
                'only mode "luminosity_density"' 'and "flux" are implemented'
            )
