import os
import sys
import argparse
import warnings

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class ArepoSnapshot:
    def __init__(
        self,
        filename,
        species,
        speciesfile,
        alpha=0.0,
        beta=0.0,
        gamma=0.0,
        boxsize=1e12,
        resolution=512,
        numthreads=4,
    ):
        """
        Loads relevant data for conversion from Arepo snapshot to a
        csvy-model. Requires arepo-snap-util to be installed.
        The snapshot is mapped onto a carthesian grid before further
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
        boxsize : float
            Size of the box (in cm) from which data is mapped
            to a Carthesian grid. Only usable with snapshots.
            Default: 1e12
        resolution : int
            Resolution of the Carthesian grid. Only usable
            with snapshots. Default: 512
        numthreads : int
            Number of threads with which Carthesian mapping
            is done. Default: 4
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
            filename, hdf5=True, quiet=True, lazy_load=True,
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

        self.pos = np.array(
            self.s.mapOnCartGrid(
                "pos",
                box=[boxsize, boxsize, boxsize],
                center=self.s.centerofmass(),
                res=resolution,
                numthreads=numthreads,
            )
        )
        for i in range(3):
            self.pos[i] -= self.s.centerofmass()[i]

        self.rho = np.array(
            self.s.mapOnCartGrid(
                "rho",
                box=[boxsize, boxsize, boxsize],
                center=self.s.centerofmass(),
                res=resolution,
                numthreads=numthreads,
            )
        )

        self.vel = np.array(
            self.s.mapOnCartGrid(
                "vel",
                box=[boxsize, boxsize, boxsize],
                center=self.s.centerofmass(),
                res=resolution,
                numthreads=numthreads,
            )
        )

        self.nuc_dict = {}

        for i, spec in enumerate(self.species):
            self.nuc_dict[spec] = np.array(
                self.nucMapOnCartGrid(
                    self.s,
                    spec,
                    self.spec_ind[i],
                    box=[boxsize, boxsize, boxsize],
                    res=resolution,
                    center=self.s.centerofmass(),
                    numthreads=numthreads,
                )
            )

    def nucMapOnCartGrid(
        self,
        snapshot,
        species,
        ind,
        box,
        res=512,
        numthreads=1,
        value="xnuc",
        center=False,
        saveas=False,
        use_only_cells=None,
    ):
        """
        Helper funciton to extract nuclear composition from snapshots
        """

        try:
            import pylab
            import calcGrid
        except ModuleNotFoundError:
            raise ImportError(
                "Please make sure you have arepo-snap-util installed if you want to directly import Arepo snapshots."
            )
        if type(center) == list:
            center = pylab.array(center)
        elif type(center) != np.ndarray:
            center = snapshot.center

        if type(box) == list:
            box = pylab.array(box)
        elif type(box) != np.ndarray:
            box = np.array(
                [snapshot.boxsize, snapshot.boxsize, snapshot.boxsize]
            )

        if type(res) == list:
            res = pylab.array(res)
        elif type(res) != np.ndarray:
            res = np.array([res] * 3)

        if use_only_cells is None:
            use_only_cells = np.arange(snapshot.nparticlesall[0], dtype="int32")

        pos = snapshot.pos[use_only_cells, :].astype("float64")
        px = np.abs(pos[:, 0] - center[0])
        py = np.abs(pos[:, 1] - center[1])
        pz = np.abs(pos[:, 2] - center[2])

        (pp,) = np.where(
            (px < 0.5 * box[0]) & (py < 0.5 * box[1]) & (pz < 0.5 * box[2])
        )
        print("Selected %d of %d particles." % (pp.size, snapshot.npart))

        posdata = pos[pp]
        valdata = snapshot.data[value][use_only_cells, ind][pp].astype(
            "float64"
        )

        if valdata.ndim == 1:
            data = calcGrid.calcASlice(
                posdata,
                valdata,
                nx=res[0],
                ny=res[1],
                nz=res[2],
                boxx=box[0],
                boxy=box[1],
                boxz=box[2],
                centerx=center[0],
                centery=center[1],
                centerz=center[2],
                grid3D=True,
                numthreads=numthreads,
            )
            grid = data["grid"]
        else:
            # We are going to generate ndim 3D grids and stack them together
            # in a grid of shape (valdata.shape[1],res,res,res)
            grid = []
            for dim in range(valdata.shape[1]):
                data = calcGrid.calcASlice(
                    posdata,
                    valdata[:, dim],
                    nx=res[0],
                    ny=res[1],
                    nz=res[2],
                    boxx=box[0],
                    boxy=box[1],
                    boxz=box[2],
                    centerx=center[0],
                    centery=center[1],
                    centerz=center[2],
                    grid3D=True,
                    numthreads=numthreads,
                )
                grid.append(data["grid"])
            grid = np.stack([subgrid for subgrid in grid])
        if saveas:
            grid.tofile(saveas)

        return grid

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
        -----
        pos : list of float
            Meshgrid of positions in center of mass frames in
            Carthesian coordinates
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
        -----
        save : str
            Path under which the figure is to be saved. Default: None
        dpi : int
            Dpi of the saved figure
        **kwargs : keywords passable to matplotlib.pyplot.plot()

        Returns
        -----
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
        ax2.set_xlabel("Radial position (cm)")  # TODO astropy unit support
        ax2.set_title("Profiles along the positive axis")

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
                save, bbox_inches="tight", dpi=dpi,
            )

        return fig

    def export(self, nshells, filename, direction="pos"):
        """
        Function to export a profile as csvy file. Either the
        positive or negative direction can be exported. Does
        not overwrite existing files, saves to *(<n>).csvy
        file instead.

        Parameters
        -----
        nshells : int
            Number of shells to be exported.
        filename : str
            Name of the exported file
        direction : str
            Specifies if either the positive or negative
            direction is to be exported. Available
            options: ['pos', 'neg']. Default: pos

        Returns
        -----
        filename : str
            Name of the actual saved file
        """

        # Find a free filename
        if filename.endswith(".csvy"):
            filename = filename.replace(".csvy", "")

        if os.path.exists("%s.csvy" % filename):
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

            f.write("".join(["\n", "---\n",]))

            # WRITE DATA
            datastring = ["velocity,", "density,"]
            for spec in self.species[:-1]:
                datastring.append("%s," % spec.capitalize())
            datastring.append("%s" % self.species[-1].capitalize())
            f.write("".join(datastring))

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

            # Selcect nshells equally spaced indices
            if nshells > len(exp):
                warnings.warn(
                    "nshells was grater then available resolution. Setting nshells to resolution"
                )
                nshells = len(exp)
            inds = np.linspace(0, len(exp[0]) - 1, num=nshells, dtype=int)

            for i in inds:
                f.write("\n")
                for ii in range(len(exp) - 1):
                    f.write("%g," % exp[ii][i])
                f.write("%g" % exp[-1][i])

        return filename


class LineProfile(Profile):
    """
    Class for profiles extrected along a line, i.e. the x-axis.
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
        Creates a profile along the x-axis without any averaging

        Parameters
        -----
        inner_radius : float
            Inner radius where the profiles will be cut off. Default: None
        outer_radius : float
            Outer radius where the profiles will be cut off. Default: None
        show_plot : bool
            Specifies if a plot is to be shown after the creation of the
            profile. Default: True
        save_plot : str
            Location where the plot is being saved. Default: None
        plot_dpi : int
            Dpi of the saved plot. Default: 600

        Returns
        -----
        profile : LineProfile object

        """

        midpoint = int(np.ceil(len(self.rho) / 2))

        # Extract radialprofiles
        pos_p = np.sqrt(
            (self.pos[0, midpoint, midpoint:, midpoint]) ** 2
            + (self.pos[1, midpoint, midpoint:, midpoint]) ** 2
            + (self.pos[2, midpoint, midpoint:, midpoint]) ** 2
        )
        pos_n = np.sqrt(
            self.pos[0, midpoint, :midpoint, midpoint] ** 2
            + self.pos[1, midpoint, :midpoint, midpoint] ** 2
            + self.pos[2, midpoint, :midpoint, midpoint] ** 2
        )

        vel_p = np.sqrt(
            self.vel[0, midpoint, midpoint:, midpoint] ** 2
            + self.vel[1, midpoint, midpoint:, midpoint] ** 2
            + self.vel[2, midpoint, midpoint:, midpoint] ** 2
        )
        vel_n = np.sqrt(
            self.vel[0, midpoint, :midpoint, midpoint] ** 2
            + self.vel[1, midpoint, :midpoint, midpoint] ** 2
            + self.vel[2, midpoint, :midpoint, midpoint] ** 2
        )

        rho_p = self.rho[midpoint, midpoint:, midpoint]
        rho_n = self.rho[midpoint, :midpoint, midpoint]

        spec_p = {}
        spec_n = {}

        for spec in self.species:
            spec_p[spec] = self.xnuc[spec][midpoint, midpoint:, midpoint]
            spec_n[spec] = self.xnuc[spec][midpoint, :midpoint, midpoint]

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

        if show_plot:
            self.plot_profile(save=save_plot, dpi=plot_dpi)

        return self


class ConeProfile(Profile):
    pass


class FullProfile(Profile):
    pass
