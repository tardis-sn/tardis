""" This module provides a line identification tool that allows the user to
extract a summary of the most prominent transitions in a TARDIS simulation. """

import numpy as np
import pandas as pd
import astropy.units as u

import matplotlib.pyplot as plt
import matplotlib.ticker as tck

from tardis.util.base import (
    atomic_number2element_symbol,
    int_to_roman,
)

# Import the SDECData class from the sdec_plot module
from tardis.visualization.tools.sdec_plot import SDECData


class LineIdentifier(object):
    def __init__(self, data):
        """
        Initialize line_identifier with required data of simulation model.

        Parameters
        ----------
        data : dict of SDECData
            Dictionary to store data required for SDEC plot, for both packet
            modes i.e. real and virtual
        """
        self.data = data

    @classmethod
    def from_simulation(cls, sim):
        """
        Create an instance of line_identifier from a TARDIS simulation object.

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            TARDIS Simulation object produced by running a simulation

        Returns
        -------
        line_identifier
        """
        if sim.transport.enable_vpacket_tracking:
            return cls(SDECData.from_simulation(sim, "virtual"))
        else:
            return cls(SDECData.from_simulation(sim, "real"))

    def _reset_cache(self):
        self._reset_lam_min()
        self._reset_lam_max()
        self._reset_derived_quantities()

    def _reset_lam_min(self):
        self._lam_min = None

    def _reset_lam_max(self):
        self._lam_max = None

    def _reset_derived_quantities(self):
        self._line_mask = None
        self._lam_in = None
        self._lines_info_unique = None
        self._lines_count = None
        self._lines_ids = None
        self._lines_ids_unique = None

    @property
    def lam_min(self):
        if self._lam_min is None:
            raise ValueError("lam_min not set")
        return self._lam_min

    @lam_min.setter
    def lam_min(self, val):
        self._reset_derived_quantities()
        try:
            self._lam_min = val.to(u.AA)
        except AttributeError:
            self._lam_min = val * u.AA

    @property
    def lam_max(self):
        if self._lam_max is None:
            raise ValueError("lam_max not set")
        return self._lam_max

    @lam_max.setter
    def lam_max(self, val):
        self._reset_derived_quantities()
        try:
            self._lam_max = val.to(u.AA)
        except AttributeError:
            self._lam_max = val * u.AA

    @property
    def lam_in(self):
        if self._lam_in is None:
            self._lam_in = u.Quantity(
                self.data.packets_df.last_line_interaction_in_nu, u.Hz
            ).to(u.AA, equivalencies=u.spectral())
        return self._lam_in

    @property
    def line_mask(self):
        if self._line_mask is None:
            self._line_mask = (self.lam_in >= self.lam_min) * (
                self.lam_in <= self.lam_max
            )
        return self._line_mask

    @property
    def lines_ids(self):
        if self._lines_ids is None:
            ids = self.data.packets_df.last_line_interaction_in_id[
                self.line_mask
            ].values
            self._lines_ids = self.data.lines_df.iloc[ids].index
        return self._lines_ids

    @property
    def lines_ids_unique(self):
        if self._lines_ids_unique is None:
            self._lines_ids_unique = np.unique(self.lines_ids)
        return self._lines_ids_unique

    @property
    def lines_info_unique(self):
        if self._lines_info_unique is None:
            self._lines_info_unique = self.data.lines_df.loc[
                self.lines_ids_unique
            ]
        return self._lines_info_unique

    @property
    def lines_count(self):
        if self._lines_count is None:
            counts = np.bincount(self.lines_ids)
            self._lines_count = counts[counts > 0]
        return self._lines_count

    def identify(self, lam_min, lam_max):
        self.lam_min = lam_min
        self.lam_max = lam_max

    def plot_summary(
        self, nlines=None, lam_min=None, lam_max=None, output_filename=None
    ):

        fig = plt.figure()
        ax = fig.add_subplot()
        fig.subplots_adjust(left=0.2)

        ax.tick_params(
            axis="x", direction="in", top=True, bottom=True, which="both"
        )
        ax.xaxis.set_minor_locator(tck.AutoMinorLocator())

        if lam_min is None:
            self.lam_min = np.min(self.data.spectrum_wavelength).value
        else:
            self.lam_min = lam_min

        if lam_max is None:
            self.lam_max = np.max(self.data.spectrum_wavelength).value
        else:
            self.lam_max = lam_max

        _lines_count = self.lines_count[np.argsort(self.lines_count)][::-1]

        _lines_fraction = self.lines_count[np.argsort(self.lines_count)][
            ::-1
        ] / float(self.lines_count.sum())

        _lines_ids = self.lines_ids_unique[np.argsort(self.lines_count)][::-1]

        if nlines is None:
            if len(_lines_count) > 20:
                self.nlines = 20
            else:
                self.nlines = len(_lines_count)
        else:
            if len(_lines_count) > nlines:
                self.nlines = nlines
            else:
                self.nlines = len(_lines_count)

        species = []
        wavelengths = []
        labels = []

        for line_id in _lines_ids:
            chemical_symbol = atomic_number2element_symbol(
                self.lines_info_unique.loc[line_id].atomic_number
            )

            ionisation_level = int_to_roman(
                int(self.lines_info_unique.loc[line_id].ion_number) + 1
            )

            species.append(f"{chemical_symbol} {ionisation_level}")

            wavelengths.append(
                f"{self.lines_info_unique.loc[line_id].wavelength:.3f}"
            )

            labels.append(
                f"{chemical_symbol}$\,${ionisation_level} $\lambda${self.lines_info_unique.loc[line_id].wavelength:.3f}"
            )

        # Parameters for the output plot
        ax.set_title(
            f"Line transitions for ${self.lam_min.value:.1f} \leq \lambda (\mathrm{{\AA}}) \leq {self.lam_max.value:.1f}$"
        )

        ax.barh(np.arange(self.nlines), _lines_fraction[: self.nlines][::-1])

        ax.set_xlabel("Fraction of total line transitions in wavelength range")

        ax.set_yticks(np.arange(len(_lines_fraction[: self.nlines][::-1])))
        ax.set_yticklabels(labels[: self.nlines][::-1], size="medium")

        ax.annotate(
            f"{self.nlines}/{len(self.lines_count)} lines displayed\n {np.sum(_lines_count[:self.nlines])}/{len(self.lines_ids)} interacting and\n escaping r-packets displayed",
            xy=(0.95, 0.05),
            xycoords="axes fraction",
            horizontalalignment="right",
            verticalalignment="bottom",
        )

        """ Below we export the compiled line information to output .csv files,
        IF the user specifies an output filename for them. """

        # This gives output file for contributions of lines
        if output_filename != None:
            dataframe = pd.DataFrame(
                {
                    "Species": species,
                    "Wavelength(AA)": wavelengths,
                    "Total_no_transitions": _lines_count,
                    "Fraction_total_transitions": _lines_fraction,
                }
            )

            with open(output_filename, "w") as file:
                file.write(
                    f"# Line transitions in wavelength range {self.lam_min.value:.1f}-{self.lam_max.value:.1f} Angstroms \n"
                )

                dataframe.to_csv(file, sep="\t", index=False)

        # This gives collapsed contributions for species summed across all lines
        if output_filename != None:
            dataframe = pd.DataFrame(
                {
                    "Species": species,
                    "Wavelength(AA)": wavelengths,
                    "Total_no_transitions": _lines_count,
                    "Fraction_total_transitions": _lines_fraction,
                }
            )

            dataframe_collapsed = dataframe.groupby(["Species"]).sum()
            dataframe_collapsed = dataframe_collapsed.drop(
                columns=["Wavelength(AA)"]
            )

            dataframe_collapsed = dataframe_collapsed.sort_values(
                by=["Total_no_transitions"], ascending=False
            )

            with open(
                f"{output_filename[:-4]}_collapsed{output_filename[-4:]}", "w"
            ) as file:
                file.write(
                    f"# Line transitions in wavelength range {self.lam_min.value:.1f}-{self.lam_max.value:.1f} Angstroms \n"
                )

                dataframe_collapsed.to_csv(file, sep="\t", index=True)
