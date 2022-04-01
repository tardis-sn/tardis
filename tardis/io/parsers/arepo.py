import os
import sys
import argparse
import warnings

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats


class ArepoSnapshot:
    def __init__(
        self,
        filename,
        species,
        speciesfile,
        alpha=0.0,
        beta=0.0,
        gamma=0.0,
    ):
        """
        Loads relevant data for conversion from Arepo snapshot to a
        csvy-model. Requires arepo-snap-util to be installed.
        The snapshot is mapped onto a Cartesian grid before further
        processing is done.

        Parameters
        ----------
        filename : str
            Path to file to be converted.
        species : list of str
            Names of the species to be exported. Have to be the
            same as in the species-file of the Arepo simulation
        speciesfile : str
            File specifying the species used in the Arepo
            simulation.
        alpha : float
            Euler angle alpha for rotation of the desired line-
            of-sight to the x-axis. Only usable with snapshots.
            Default: 0.0
        beta : float
            Euler angle beta for rotation of the desired line-
            of-sight to the x-axis. Only usable with snapshots.
            Default: 0.0
        gamma : float
            Euler angle gamma for rotation of the desired line-
            of-sight to the x-axis. Only usable with snapshots.
            Default: 0.0
        """

        try:
            import gadget_snap
            import calcGrid
        except ModuleNotFoundError:
            raise ImportError(
                "Please make sure you have arepo-snap-util installed if you want to directly import Arepo snapshots."
            )

        self.species = species
        species_full = np.genfromtxt(speciesfile, skip_header=1, dtype=str).T[0]
        self.spec_ind = []
        for spec in self.species:
            self.spec_ind.append(np.where(species_full == spec)[0][0])

        self.spec_ind = np.array(self.spec_ind)

        self.s = gadget_snap.gadget_snapshot(
            filename,
            hdf5=True,
            quiet=True,
            lazy_load=True,
        )

        rz_yaw = np.array(
            [
                [np.cos(alpha), -np.sin(alpha), 0],
                [np.sin(alpha), np.cos(alpha), 0],
                [0, 0, 1],
            ]
        )
        ry_pitch = np.array(
            [
                [np.cos(beta), 0, np.sin(beta)],
                [0, 1, 0],
                [-np.sin(beta), 0, np.cos(beta)],
            ]
        )
        rx_roll = np.array(
            [
                [1, 0, 0],
                [0, np.cos(gamma), -np.sin(gamma)],
                [0, np.sin(gamma), np.cos(gamma)],
            ]
        )
        # R = RzRyRx
        rotmat = np.dot(rz_yaw, np.dot(ry_pitch, rx_roll))

        self.s.rotateto(rotmat[0], dir2=rotmat[1], dir3=rotmat[2])

        self.time = self.s.time
        self.pos = np.array(self.s.data["pos"][: self.s.nparticlesall[0]])
        self.pos = self.pos.T
        # Update position to CoM frame
        for i in range(3):
            self.pos[i] -= self.s.centerofmass()[i]
        self.rho = np.array(self.s.data["rho"])
        self.vel = np.array(self.s.data["vel"][: self.s.nparticlesall[0]])
        self.vel = self.vel.T
        self.nuc_dict = {}

        for i, spec in enumerate(self.species):
            self.nuc_dict[spec] = np.array(
                self.s.data["xnuc"][:, self.spec_ind[i]]
            )

    def get_grids(self):
        """
        Returns all relevant data to create Profile objects
        """
        return self.pos, self.vel, self.rho, self.nuc_dict, self.time


class Profile:
    """
    Parent class of all Profiles. Contains general function,
    e.g. for plotting and export.
    """

    def __init__(self, pos, vel, rho, xnuc, time):
        """
        Parameters
        ----------
        pos : list of float
            Meshgrid of positions in center of mass frames in
            Cartesian coordinates
        vel : list of float
            Meshgrid of velocities/ velocity vectors
        rho : list of float
            Meshgrid of density
        xnuc : dict
            Dictonary containing all the nuclear fraction
            meshgrids of the relevant species.
        time : float
            Time of the data

        """

        self.pos = pos
        self.vel = vel
        self.rho = rho
        self.xnuc = xnuc
        self.time = time

        self.species = list(self.xnuc.keys())

        # Empty values to be filled with the create_profile function
        self.pos_prof_p = None
        self.pos_prof_n = None

        self.vel_prof_p = None
        self.vel_prof_n = None

        self.rho_prof_p = None
        self.rho_prof_n = None

        self.xnuc_prof_p = {}
        self.xnuc_prof_n = {}

    def plot_profile(self, save=None, dpi=600, **kwargs):
        """
        Plots profile, both in the positive and negative direction.

        Parameters
        ----------
        save : str
            Path under which the figure is to be saved. Default: None
        dpi : int
            Dpi of the saved figure
        **kwargs : keywords passable to matplotlib.pyplot.plot()

        Returns
        -------
        fig : matplotlib figure object
        """

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=[9.8, 9.6])

        # Positive direction plots
        ax1.plot(
            self.pos_prof_p,
            self.rho_prof_p / max(self.rho_prof_p),
            label="Density",
            **kwargs,
        )
        ax1.plot(
            self.pos_prof_p,
            self.vel_prof_p / max(self.vel_prof_p),
            label="Velocity",
            **kwargs,
        )
        for spec in self.species:
            ax1.plot(
                self.pos_prof_p,
                self.xnuc_prof_p[spec],
                label=spec.capitalize(),
                **kwargs,
            )

        ax1.grid()
        ax1.set_ylabel("Profile (arb. unit)")
        ax1.set_title("Profiles along the positive axis")

        # Positive direction plots
        ax2.plot(
            self.pos_prof_n,
            self.rho_prof_n / max(self.rho_prof_n),
            label="Density",
            **kwargs,
        )
        ax2.plot(
            self.pos_prof_n,
            self.vel_prof_n / max(self.vel_prof_n),
            label="Velocity",
            **kwargs,
        )
        for spec in self.species:
            ax2.plot(
                self.pos_prof_n,
                self.xnuc_prof_n[spec],
                label=spec.capitalize(),
                **kwargs,
            )

        ax2.grid()
        ax2.set_ylabel("Profile (arb. unit)")
        ax2.set_xlabel("Radial position (cm)")
        ax2.set_title("Profiles along the negative axis")

        # Some styling
        fig.tight_layout()

        handles, labels = ax1.get_legend_handles_labels()
        lgd = ax1.legend(
            handles,
            labels,
            loc="upper left",
            bbox_to_anchor=(1.05, 1.05),
            title="Time = {:.2f} s".format(self.time),
        )
        if save is not None:
            plt.savefig(
                save,
                bbox_inches="tight",
                dpi=dpi,
            )

        return fig

    def rebin(self, nshells, statistic="mean"):
        """
        Rebins the data to nshells. Uses the scipy.stats.binned_statistic
        to bin the data. The standard deviation of each bin can be obtained
        by passing the statistics="std" keyword.

        Parameters
        ----------
        nshells : int
            Number of bins of new data.
        statistic : str
            Scipy keyword for scipy.stats.binned_statistic. Default: mean

        Returns
        -------
        self : Profile object

        """

        self.vel_prof_p, bins_p = stats.binned_statistic(
            self.pos_prof_p,
            self.vel_prof_p,
            statistic=statistic,
            bins=nshells,
        )[:2]
        self.vel_prof_n, bins_n = stats.binned_statistic(
            self.pos_prof_n,
            self.vel_prof_n,
            statistic=statistic,
            bins=nshells,
        )[:2]

        self.rho_prof_p = stats.binned_statistic(
            self.pos_prof_p,
            self.rho_prof_p,
            statistic=statistic,
            bins=nshells,
        )[0]
        self.rho_prof_n = stats.binned_statistic(
            self.pos_prof_n,
            self.rho_prof_n,
            statistic=statistic,
            bins=nshells,
        )[0]

        for spec in self.species:
            self.xnuc_prof_p[spec] = stats.binned_statistic(
                self.pos_prof_p,
                self.xnuc_prof_p[spec],
                statistic=statistic,
                bins=nshells,
            )[0]
            self.xnuc_prof_n[spec] = stats.binned_statistic(
                self.pos_prof_n,
                self.xnuc_prof_n[spec],
                statistic=statistic,
                bins=nshells,
            )[0]

        self.pos_prof_p = np.array(
            [(bins_p[i] + bins_p[i + 1]) / 2 for i in range(len(bins_p) - 1)]
        )
        self.pos_prof_n = np.array(
            [(bins_n[i] + bins_n[i + 1]) / 2 for i in range(len(bins_n) - 1)]
        )

        return self

    def export(
        self,
        nshells,
        filename,
        direction="pos",
        statistic="mean",
        overwrite=False,
    ):
        """
        Function to export a profile as csvy file. Either the
        positive or negative direction can be exported. By default
        does not overwrite existing files, saves to <filename>_<number>.csvy
        file instead.

        Parameters
        ----------
        nshells : int
            Number of shells to be exported.
        filename : str
            Name of the exported file
        direction : str
            Specifies if either the positive or negative
            direction is to be exported. Available
            options: ['pos', 'neg']. Default: pos
        statistic : str
            Scipy keyword for scipy.stats.binned_statistic. If
            statistic=None, data is not rebinned. Default: "mean"
        overwrite: bool
            If true, will overwrite if a file of the same name exists.
            By default False.

        Returns
        -------
        filename : str
            Name of the actual saved file
        """

        # Find a free filename
        if filename.endswith(".csvy"):
            filename = filename.replace(".csvy", "")

        if os.path.exists("%s.csvy" % filename) and not overwrite:
            i = 0
            while os.path.exists("%s_%s.csvy" % (filename, i)):
                i += 1
            filename = "%s_%s.csvy" % (filename, i)
        else:
            filename = "%s.csvy" % filename

        with open(filename, "w") as f:
            # WRITE HEADER
            f.write(
                "".join(
                    [
                        "---\n",
                        "name: csvy_full\n",
                        "model_density_time_0: {:g} day\n".format(
                            self.time / (3600 * 24)
                        ),  # TODO astropy units
                        "model_isotope_time_0: {:g} day\n".format(
                            self.time / (3600 / 24)
                        ),  # TODO astropy units
                        "description: Config file for TARDIS from Arepo snapshot.\n",
                        "tardis_model_config_version: v1.0\n",
                        "datatype:\n",
                        "  fields:\n",
                        "    -  name: velocity\n",
                        "       unit: cm/s\n",
                        "       desc: velocities of shell outer bounderies.\n",
                        "    -  name: density\n",
                        "       unit: g/cm^3\n",
                        "       desc: density of shell.\n",
                    ]
                )
            )

            for spec in self.species:
                f.write(
                    "".join(
                        [
                            "    -  name: %s\n" % spec.capitalize(),
                            "       desc: fractional %s abundance.\n"
                            % spec.capitalize(),
                        ]
                    )
                )

            f.write(
                "".join(
                    [
                        "\n",
                        "---\n",
                    ]
                )
            )

            # WRITE DATA
            datastring = ["velocity,", "density,"]
            for spec in self.species[:-1]:
                datastring.append("%s," % spec.capitalize())
            datastring.append("%s" % self.species[-1].capitalize())
            f.write("".join(datastring))

            # Rebin data to nshells
            if statistic is not None:
                self.rebin(nshells, statistic=statistic)

            if direction == "pos":
                exp = [
                    self.vel_prof_p,
                    self.rho_prof_p,
                ]
                for spec in self.xnuc_prof_p:
                    exp.append(self.xnuc_prof_p[spec])
            elif direction == "neg":
                exp = [
                    self.vel_prof_n,
                    self.rho_prof_n,
                ]
                for spec in self.xnuc_prof_n:
                    exp.append(self.xnuc_prof_n[spec])
            else:
                raise ValueError("Unrecognized option for keyword 'direction'")

            inds = np.linspace(0, len(exp[0]) - 1, num=nshells, dtype=int)

            for i in inds:
                f.write("\n")
                for ii in range(len(exp) - 1):
                    f.write("%g," % exp[ii][i])
                f.write("%g" % exp[-1][i])

        return filename

    def get_profiles(self):
        """Returns all profiles for manual post_processing etc."""
        return (
            self.pos_prof_p,
            self.pos_prof_n,
            self.vel_prof_p,
            self.vel_prof_n,
            self.rho_prof_p,
            self.rho_prof_n,
            self.xnuc_prof_p,
            self.xnuc_prof_n,
        )


class ConeProfile(Profile):
    """
    Class for profiles extracted inside a cone around the x-axis.
    Extends Profile.
    """

    def create_profile(
        self,
        opening_angle=20.0,
        inner_radius=None,
        outer_radius=None,
        show_plot=True,
        save_plot=None,
        plot_dpi=600,
    ):
        """
        Creates a profile along the x-axis without any averaging

        Parameters
        ----------
        opening_angle : float
            Opening angle (in degrees) of the cone from which the
            data is extracted. Refers to the total opening angle, not
            the angle with respect to the x axis. Default: 20.0
        inner_radius : float
            Inner radius where the profiles will be cut off. Default: None
        outer_radius : float
            Outer radius where the profiles will be cut off. Default: None

        Returns
        -------
        profile : ConeProfile object

        """

        # Convert Cartesian coordinates into cylindrical coordinates
        # P(x,y,z) -> P(x,r,theta)
        cyl = np.array(
            [
                self.pos[0],
                np.sqrt(self.pos[1] ** 2 + self.pos[2] ** 2),
                np.arctan(self.pos[2] / self.pos[1]),
            ]
        )

        # Get maximum allowed r of points to still be in cone
        dist = np.tan(opening_angle / 2) * np.abs(cyl[0])

        # Create masks
        cmask_p = np.logical_and(cyl[0] > 0, cyl[1] <= dist)
        cmask_n = np.logical_and(cyl[0] < 0, cyl[1] <= dist)

        # Apply mask to data
        pos_p = np.sqrt(
            (self.pos[0][cmask_p]) ** 2
            + (self.pos[1][cmask_p]) ** 2
            + (self.pos[2][cmask_p]) ** 2
        )
        pos_n = np.sqrt(
            self.pos[0][cmask_n] ** 2
            + self.pos[1][cmask_n] ** 2
            + self.pos[2][cmask_n] ** 2
        )

        vel_p = np.sqrt(
            self.vel[0][cmask_p] ** 2
            + self.vel[1][cmask_p] ** 2
            + self.vel[2][cmask_p] ** 2
        )
        vel_n = np.sqrt(
            self.vel[0][cmask_n] ** 2
            + self.vel[1][cmask_n] ** 2
            + self.vel[2][cmask_n] ** 2
        )

        rho_p = self.rho[cmask_p]
        rho_n = self.rho[cmask_n]

        spec_p = {}
        spec_n = {}

        for spec in self.species:
            spec_p[spec] = self.xnuc[spec][cmask_p]
            spec_n[spec] = self.xnuc[spec][cmask_n]

        self.pos_prof_p = np.sort(pos_p)
        self.pos_prof_n = np.sort(pos_n)

        if outer_radius is None:
            maxradius_p = max(self.pos_prof_p)
            maxradius_n = max(self.pos_prof_n)
        else:
            maxradius_p = outer_radius
            maxradius_n = outer_radius

        if inner_radius is None:
            minradius_p = min(self.pos_prof_p)
            minradius_n = min(self.pos_prof_n)
        else:
            minradius_p = inner_radius
            minradius_n = inner_radius

        mask_p = np.logical_and(
            self.pos_prof_p >= minradius_p, self.pos_prof_p <= maxradius_p
        )
        mask_n = np.logical_and(
            self.pos_prof_n >= minradius_n, self.pos_prof_n <= maxradius_n
        )

        if not mask_p.any() or not mask_n.any():
            raise ValueError("No points left between inner and outer radius.")

        self.rho_prof_p = np.array(
            [x for _, x in sorted(zip(pos_p, rho_p), key=lambda pair: pair[0])]
        )[mask_p]
        self.rho_prof_n = np.array(
            [x for _, x in sorted(zip(pos_n, rho_n), key=lambda pair: pair[0])]
        )[mask_n]

        self.vel_prof_p = np.array(
            [x for _, x in sorted(zip(pos_p, vel_p), key=lambda pair: pair[0])]
        )[mask_p]
        self.vel_prof_n = np.array(
            [x for _, x in sorted(zip(pos_n, vel_n), key=lambda pair: pair[0])]
        )[mask_n]

        for spec in self.species:
            self.xnuc_prof_p[spec] = np.array(
                [
                    x
                    for _, x in sorted(
                        zip(pos_p, spec_p[spec]), key=lambda pair: pair[0]
                    )
                ]
            )[mask_p]
            self.xnuc_prof_n[spec] = np.array(
                [
                    x
                    for _, x in sorted(
                        zip(pos_n, spec_n[spec]), key=lambda pair: pair[0]
                    )
                ]
            )[mask_n]

        self.pos_prof_p = self.pos_prof_p[mask_p]
        self.pos_prof_n = self.pos_prof_n[mask_n]

        return self


class FullProfile(Profile):
    """
    Class for profiles extracted from the full snapshot,
    i.e. angle averaged profiles.
    Extends Profile.
    """

    def create_profile(
        self,
        inner_radius=None,
        outer_radius=None,
        show_plot=True,
        save_plot=None,
        plot_dpi=600,
    ):
        """
        Creates a profile from the full snapshot. Positive and negative
        direction are identical.

        Parameters
        ----------
        inner_radius : float
            Inner radius where the profiles will be cut off. Default: None
        outer_radius : float
            Outer radius where the profiles will be cut off. Default: None

        Returns
        -------
        profile : FullProfile object

        """

        pos_p = np.sqrt(
            (self.pos[0]) ** 2 + (self.pos[1]) ** 2 + (self.pos[2]) ** 2
        ).flatten()
        pos_n = np.sqrt(
            self.pos[0] ** 2 + self.pos[1] ** 2 + self.pos[2] ** 2
        ).flatten()

        vel_p = np.sqrt(
            self.vel[0] ** 2 + self.vel[1] ** 2 + self.vel[2] ** 2
        ).flatten()
        vel_n = np.sqrt(
            self.vel[0] ** 2 + self.vel[1] ** 2 + self.vel[2] ** 2
        ).flatten()

        rho_p = self.rho.flatten()
        rho_n = self.rho.flatten()

        spec_p = {}
        spec_n = {}

        for spec in self.species:
            spec_p[spec] = self.xnuc[spec].flatten()
            spec_n[spec] = self.xnuc[spec].flatten()

        self.pos_prof_p = np.sort(pos_p)
        self.pos_prof_n = np.sort(pos_n)

        if outer_radius is None:
            maxradius_p = max(self.pos_prof_p)
            maxradius_n = max(self.pos_prof_n)
        else:
            maxradius_p = outer_radius
            maxradius_n = outer_radius

        if inner_radius is None:
            minradius_p = min(self.pos_prof_p)
            minradius_n = min(self.pos_prof_n)
        else:
            minradius_p = inner_radius
            minradius_n = inner_radius

        mask_p = np.logical_and(
            self.pos_prof_p >= minradius_p, self.pos_prof_p <= maxradius_p
        )
        mask_n = np.logical_and(
            self.pos_prof_n >= minradius_n, self.pos_prof_n <= maxradius_n
        )

        if not mask_p.any() or not mask_n.any():
            raise ValueError("No points left between inner and outer radius.")

        self.rho_prof_p = np.array(
            [x for _, x in sorted(zip(pos_p, rho_p), key=lambda pair: pair[0])]
        )[mask_p]
        self.rho_prof_n = np.array(
            [x for _, x in sorted(zip(pos_n, rho_n), key=lambda pair: pair[0])]
        )[mask_n]

        self.vel_prof_p = np.array(
            [x for _, x in sorted(zip(pos_p, vel_p), key=lambda pair: pair[0])]
        )[mask_p]
        self.vel_prof_n = np.array(
            [x for _, x in sorted(zip(pos_n, vel_n), key=lambda pair: pair[0])]
        )[mask_n]

        for spec in self.species:
            self.xnuc_prof_p[spec] = np.array(
                [
                    x
                    for _, x in sorted(
                        zip(pos_p, spec_p[spec]), key=lambda pair: pair[0]
                    )
                ]
            )[mask_p]
            self.xnuc_prof_n[spec] = np.array(
                [
                    x
                    for _, x in sorted(
                        zip(pos_n, spec_n[spec]), key=lambda pair: pair[0]
                    )
                ]
            )[mask_n]

        self.pos_prof_p = self.pos_prof_p[mask_p]
        self.pos_prof_n = self.pos_prof_n[mask_n]

        return self


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "snapshot",
        help="Snapshot file for which to create velocity profile plot",
    )
    parser.add_argument(
        "save",
        help="Filename of exported .csvy file",
    )
    parser.add_argument(
        "-a",
        "--alpha",
        help="Euler angle alpha for rotation of desired direction to x-axis. Default: 0",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "-b",
        "--beta",
        help="Euler angle beta for rotation of desired direction to x-axis. Default: 0",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "-g",
        "--gamma",
        help="Euler angle gamma for rotation of desired direction to x-axis. Default: 0",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "-o",
        "--opening_angle",
        help="Opening angle of the cone from which profile is extracted. Default 20.0",
        type=float,
        default=20.0,
    )
    parser.add_argument(
        "-n",
        "--nshells",
        help="Number of shells to create. Default: 10",
        type=int,
        default=10,
    )
    parser.add_argument(
        "-x",
        "--boxsize",
        help="Size of the box (in cm) from which data is extracted. Default: 1e12",
        type=float,
        default=1e12,
    )
    parser.add_argument(
        "-e",
        "--elements",
        help="List of species to be included. Default: ni56",
        default="ni56",
        nargs="+",
    )
    parser.add_argument(
        "--eosspecies",
        help="Species file including all the species used in the production of the composition file. Default: species55.txt",
        default="species55.txt",
    )
    parser.add_argument(
        "--outer_radius",
        help="Outer radius to which to build profile.",
        type=float,
    )
    parser.add_argument(
        "--inner_radius",
        help="Inner radius to which to build profile.",
        type=float,
    )
    parser.add_argument(
        "--profile",
        help="How to build profile. Available options: [cone, full]. Default: cone",
        default="cone",
        choices=["cone", "full"],
    )
    parser.add_argument("--plot", help="File name of saved plot.")
    parser.add_argument(
        "--dpi", help="Dpi of saved plot. Default: 600", type=int, default=600
    )

    args = parser.parse_args()

    snapshot = ArepoSnapshot(
        args.snapshot,
        args.elements,
        args.eosspecies,
        alpha=args.alpha,
        beta=args.beta,
        gamma=args.gamma,
        boxsize=args.boxsize,
        resolution=args.resolution,
        numthreads=args.numthreads,
    )

    pos, vel, rho, xnuc, time = snapshot.get_grids()

    if args.profile == "cone":
        profile = ConeProfile(pos, vel, rho, xnuc, time)
    elif args.profile == "full":
        profile = FullProfile(pos, vel, rho, xnuc, time)

    if args.profile == "cone":
        profile.create_profile(
            opening_angle=args.opening_angle,
            inner_radius=args.inner_radius,
            outer_radius=args.outer_radius,
        )
    else:
        profile.create_profile(
            inner_radius=args.inner_radius,
            outer_radius=args.outer_radius,
        )

    profile.export(args.nshells, args.save)

    if args.plot:
        profile.plot_profile(save=args.plot, dpi=args.dpi)
