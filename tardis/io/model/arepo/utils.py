"""Utility functions for reprojecting 3D Arepo data to 1D profiles."""

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from scipy import stats

from tardis.configuration.sorting_globals import SORTING_ALGORITHM
from tardis.io.model.arepo.data import ArepoData

if TYPE_CHECKING:
    from matplotlib.figure import Figure


def create_cone_profile(
    data: ArepoData,
    opening_angle: float = 20.0,
    inner_radius: float | None = None,
    outer_radius: float | None = None,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    dict,
    dict,
]:
    """
    Create 1D profiles from Arepo data using a cone selection.

    Extracts data from within a cone around the x-axis, creating separate
    profiles for the positive and negative directions.

    Parameters
    ----------
    data : ArepoData
        The Arepo snapshot data.
    opening_angle : float, optional
        Total opening angle of the cone in degrees. Default: 20.0
    inner_radius : float, optional
        Inner radius cutoff in cm. Default: None
    outer_radius : float, optional
        Outer radius cutoff in cm. Default: None

    Returns
    -------
    pos_prof_p : numpy.ndarray
        Position profile in positive direction (cm).
    pos_prof_n : numpy.ndarray
        Position profile in negative direction (cm).
    vel_prof_p : numpy.ndarray
        Velocity profile in positive direction (cm/s).
    vel_prof_n : numpy.ndarray
        Velocity profile in negative direction (cm/s).
    rho_prof_p : numpy.ndarray
        Density profile in positive direction (g/cm^3).
    rho_prof_n : numpy.ndarray
        Density profile in negative direction (g/cm^3).
    mass_prof_p : numpy.ndarray
        Mass profile in positive direction (g).
    mass_prof_n : numpy.ndarray
        Mass profile in negative direction (g).
    xnuc_prof_p : dict
        Nuclear fraction profiles in positive direction.
    xnuc_prof_n : dict
        Nuclear fraction profiles in negative direction.

    Raises
    ------
    ValueError
        If no points remain after applying radius cuts.
    """
    # Convert Cartesian to cylindrical coordinates
    cyl = np.array(
        [
            data.position[0],
            np.sqrt(data.position[1] ** 2 + data.position[2] ** 2),
            np.arctan(data.position[2] / data.position[1]),
        ]
    )

    # Maximum allowed r for points in cone
    dist = np.tan(np.radians(opening_angle) / 2) * np.abs(cyl[0])

    # Create masks for positive/negative directions
    cmask_p = np.logical_and(cyl[0] > 0, cyl[1] <= dist)
    cmask_n = np.logical_and(cyl[0] < 0, cyl[1] <= dist)

    # Apply masks to get radial distances
    pos_p = np.sqrt(
        data.position[0][cmask_p] ** 2
        + data.position[1][cmask_p] ** 2
        + data.position[2][cmask_p] ** 2
    )
    pos_n = np.sqrt(
        data.position[0][cmask_n] ** 2
        + data.position[1][cmask_n] ** 2
        + data.position[2][cmask_n] ** 2
    )

    # Velocities
    vel_p = np.sqrt(
        data.velocities[0][cmask_p] ** 2
        + data.velocities[1][cmask_p] ** 2
        + data.velocities[2][cmask_p] ** 2
    )
    vel_n = np.sqrt(
        data.velocities[0][cmask_n] ** 2
        + data.velocities[1][cmask_n] ** 2
        + data.velocities[2][cmask_n] ** 2
    )

    mass_p = data.mass[cmask_p]
    mass_n = data.mass[cmask_n]

    rho_p = data.densities[cmask_p]
    rho_n = data.densities[cmask_n]

    # Nuclear fractions
    spec_p = {}
    spec_n = {}
    for spec in data.species:
        spec_p[spec] = data.isotope_dict[spec][cmask_p]
        spec_n[spec] = data.isotope_dict[spec][cmask_n]

    # Sort by position
    pos_prof_p = np.sort(pos_p, kind=SORTING_ALGORITHM)
    pos_prof_n = np.sort(pos_n, kind=SORTING_ALGORITHM)

    # Apply radius cuts
    maxradius_p = max(pos_prof_p) if outer_radius is None else outer_radius
    maxradius_n = max(pos_prof_n) if outer_radius is None else outer_radius
    minradius_p = min(pos_prof_p) if inner_radius is None else inner_radius
    minradius_n = min(pos_prof_n) if inner_radius is None else inner_radius

    mask_p = np.logical_and(
        pos_prof_p >= minradius_p, pos_prof_p <= maxradius_p
    )
    mask_n = np.logical_and(
        pos_prof_n >= minradius_n, pos_prof_n <= maxradius_n
    )

    if not mask_p.any() or not mask_n.any():
        raise ValueError("No points left between inner and outer radius.")

    # Sort all quantities by position
    mass_prof_p = np.array(
        [x for _, x in sorted(zip(pos_p, mass_p), key=lambda pair: pair[0])]
    )[mask_p]
    mass_prof_n = np.array(
        [x for _, x in sorted(zip(pos_n, mass_n), key=lambda pair: pair[0])]
    )[mask_n]

    rho_prof_p = np.array(
        [x for _, x in sorted(zip(pos_p, rho_p), key=lambda pair: pair[0])]
    )[mask_p]
    rho_prof_n = np.array(
        [x for _, x in sorted(zip(pos_n, rho_n), key=lambda pair: pair[0])]
    )[mask_n]

    vel_prof_p = np.array(
        [x for _, x in sorted(zip(pos_p, vel_p), key=lambda pair: pair[0])]
    )[mask_p]
    vel_prof_n = np.array(
        [x for _, x in sorted(zip(pos_n, vel_n), key=lambda pair: pair[0])]
    )[mask_n]

    xnuc_prof_p = {}
    xnuc_prof_n = {}
    for spec in data.species:
        xnuc_prof_p[spec] = np.array(
            [
                x
                for _, x in sorted(
                    zip(pos_p, spec_p[spec]), key=lambda pair: pair[0]
                )
            ]
        )[mask_p]
        xnuc_prof_n[spec] = np.array(
            [
                x
                for _, x in sorted(
                    zip(pos_n, spec_n[spec]), key=lambda pair: pair[0]
                )
            ]
        )[mask_n]

    pos_prof_p = pos_prof_p[mask_p]
    pos_prof_n = pos_prof_n[mask_n]

    return (
        pos_prof_p,
        pos_prof_n,
        vel_prof_p,
        vel_prof_n,
        rho_prof_p,
        rho_prof_n,
        mass_prof_p,
        mass_prof_n,
        xnuc_prof_p,
        xnuc_prof_n,
    )


def create_full_profile(
    data: ArepoData,
    inner_radius: float | None = None,
    outer_radius: float | None = None,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    dict,
    dict,
]:
    """
    Create 1D profiles from full Arepo snapshot (angle-averaged).

    Creates angle-averaged profiles from all data. Positive and negative
    direction profiles are identical in this case.

    Parameters
    ----------
    data : ArepoData
        The Arepo snapshot data.
    inner_radius : float, optional
        Inner radius cutoff in cm. Default: None
    outer_radius : float, optional
        Outer radius cutoff in cm. Default: None

    Returns
    -------
    pos_prof_p : numpy.ndarray
        Position profile in positive direction (cm).
    pos_prof_n : numpy.ndarray
        Position profile in negative direction (cm).
    vel_prof_p : numpy.ndarray
        Velocity profile in positive direction (cm/s).
    vel_prof_n : numpy.ndarray
        Velocity profile in negative direction (cm/s).
    rho_prof_p : numpy.ndarray
        Density profile in positive direction (g/cm^3).
    rho_prof_n : numpy.ndarray
        Density profile in negative direction (g/cm^3).
    mass_prof_p : numpy.ndarray
        Mass profile in positive direction (g).
    mass_prof_n : numpy.ndarray
        Mass profile in negative direction (g).
    xnuc_prof_p : dict
        Nuclear fraction profiles in positive direction.
    xnuc_prof_n : dict
        Nuclear fraction profiles in negative direction.

    Raises
    ------
    ValueError
        If no points remain after applying radius cuts.
    """
    # Compute radial distances
    pos_p = np.sqrt(
        data.position[0] ** 2 + data.position[1] ** 2 + data.position[2] ** 2
    ).flatten()
    pos_n = pos_p.copy()

    vel_p = np.sqrt(
        data.velocities[0] ** 2
        + data.velocities[1] ** 2
        + data.velocities[2] ** 2
    ).flatten()
    vel_n = vel_p.copy()

    mass_p = data.mass.flatten()
    mass_n = mass_p.copy()

    rho_p = data.densities.flatten()
    rho_n = rho_p.copy()

    spec_p = {}
    spec_n = {}
    for spec in data.species:
        spec_p[spec] = data.isotope_dict[spec].flatten()
        spec_n[spec] = spec_p[spec].copy()

    # Sort by position
    pos_prof_p = np.sort(pos_p, kind=SORTING_ALGORITHM)
    pos_prof_n = np.sort(pos_n, kind=SORTING_ALGORITHM)

    # Apply radius cuts
    maxradius_p = max(pos_prof_p) if outer_radius is None else outer_radius
    maxradius_n = max(pos_prof_n) if outer_radius is None else outer_radius
    minradius_p = min(pos_prof_p) if inner_radius is None else inner_radius
    minradius_n = min(pos_prof_n) if inner_radius is None else inner_radius

    mask_p = np.logical_and(
        pos_prof_p >= minradius_p, pos_prof_p <= maxradius_p
    )
    mask_n = np.logical_and(
        pos_prof_n >= minradius_n, pos_prof_n <= maxradius_n
    )

    if not mask_p.any() or not mask_n.any():
        raise ValueError("No points left between inner and outer radius.")

    # Sort all quantities by position
    mass_prof_p = np.array(
        [x for _, x in sorted(zip(pos_p, mass_p), key=lambda pair: pair[0])]
    )[mask_p]
    mass_prof_n = np.array(
        [x for _, x in sorted(zip(pos_n, mass_n), key=lambda pair: pair[0])]
    )[mask_n]

    rho_prof_p = np.array(
        [x for _, x in sorted(zip(pos_p, rho_p), key=lambda pair: pair[0])]
    )[mask_p]
    rho_prof_n = np.array(
        [x for _, x in sorted(zip(pos_n, rho_n), key=lambda pair: pair[0])]
    )[mask_n]

    vel_prof_p = np.array(
        [x for _, x in sorted(zip(pos_p, vel_p), key=lambda pair: pair[0])]
    )[mask_p]
    vel_prof_n = np.array(
        [x for _, x in sorted(zip(pos_n, vel_n), key=lambda pair: pair[0])]
    )[mask_n]

    xnuc_prof_p = {}
    xnuc_prof_n = {}
    for spec in data.species:
        xnuc_prof_p[spec] = np.array(
            [
                x
                for _, x in sorted(
                    zip(pos_p, spec_p[spec]), key=lambda pair: pair[0]
                )
            ]
        )[mask_p]
        xnuc_prof_n[spec] = np.array(
            [
                x
                for _, x in sorted(
                    zip(pos_n, spec_n[spec]), key=lambda pair: pair[0]
                )
            ]
        )[mask_n]

    pos_prof_p = pos_prof_p[mask_p]
    pos_prof_n = pos_prof_n[mask_n]

    return (
        pos_prof_p,
        pos_prof_n,
        vel_prof_p,
        vel_prof_n,
        rho_prof_p,
        rho_prof_n,
        mass_prof_p,
        mass_prof_n,
        xnuc_prof_p,
        xnuc_prof_n,
    )


def rebin_profile(
    pos_prof_p: np.ndarray,
    pos_prof_n: np.ndarray,
    vel_prof_p: np.ndarray,
    vel_prof_n: np.ndarray,
    rho_prof_p: np.ndarray,
    rho_prof_n: np.ndarray,
    mass_prof_p: np.ndarray,
    mass_prof_n: np.ndarray,
    xnuc_prof_p: dict,
    xnuc_prof_n: dict,
    nshells: int,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    dict,
    dict,
]:
    """
    Rebin profile data to specified number of shells.

    Uses scipy.stats.binned_statistic to bin the data into nshells bins.

    Parameters
    ----------
    pos_prof_p : numpy.ndarray
        Position profile in positive direction.
    pos_prof_n : numpy.ndarray
        Position profile in negative direction.
    vel_prof_p : numpy.ndarray
        Velocity profile in positive direction.
    vel_prof_n : numpy.ndarray
        Velocity profile in negative direction.
    rho_prof_p : numpy.ndarray
        Density profile in positive direction.
    rho_prof_n : numpy.ndarray
        Density profile in negative direction.
    mass_prof_p : numpy.ndarray
        Mass profile in positive direction.
    mass_prof_n : numpy.ndarray
        Mass profile in negative direction.
    xnuc_prof_p : dict
        Nuclear fraction profiles in positive direction.
    xnuc_prof_n : dict
        Nuclear fraction profiles in negative direction.
    nshells : int
        Number of bins for rebinned data.

    Returns
    -------
    pos_prof_p : numpy.ndarray
        Rebinned position profile in positive direction.
    pos_prof_n : numpy.ndarray
        Rebinned position profile in negative direction.
    vel_prof_p : numpy.ndarray
        Rebinned velocity profile in positive direction.
    vel_prof_n : numpy.ndarray
        Rebinned velocity profile in negative direction.
    rho_prof_p : numpy.ndarray
        Rebinned density profile in positive direction.
    rho_prof_n : numpy.ndarray
        Rebinned density profile in negative direction.
    mass_prof_p : numpy.ndarray
        Rebinned mass profile in positive direction.
    mass_prof_n : numpy.ndarray
        Rebinned mass profile in negative direction.
    xnuc_prof_p : dict
        Rebinned nuclear fraction profiles in positive direction.
    xnuc_prof_n : dict
        Rebinned nuclear fraction profiles in negative direction.
    """
    species = list(xnuc_prof_p.keys())

    # Bin velocities (mass-weighted mean)
    vel_prof_p_new, bins_p = stats.binned_statistic(
        pos_prof_p,
        vel_prof_p * mass_prof_p,
        statistic="mean",
        bins=nshells,
    )[:2]
    vel_prof_p_new /= stats.binned_statistic(
        pos_prof_p, mass_prof_p, statistic="mean", bins=nshells
    )[0]

    vel_prof_n_new, bins_n = stats.binned_statistic(
        pos_prof_n,
        vel_prof_n * mass_prof_n,
        statistic="mean",
        bins=nshells,
    )[:2]
    vel_prof_n_new /= stats.binned_statistic(
        pos_prof_n, mass_prof_n, statistic="mean", bins=nshells
    )[0]

    # Bin nuclear fractions (mass-weighted mean)
    xnuc_prof_p_new = {}
    xnuc_prof_n_new = {}
    for spec in species:
        xnuc_prof_p_new[spec] = (
            stats.binned_statistic(
                pos_prof_p,
                xnuc_prof_p[spec] * mass_prof_p,
                statistic="mean",
                bins=nshells,
            )[0]
            / stats.binned_statistic(
                pos_prof_p, mass_prof_p, statistic="mean", bins=nshells
            )[0]
        )

        xnuc_prof_n_new[spec] = (
            stats.binned_statistic(
                pos_prof_n,
                xnuc_prof_n[spec] * mass_prof_n,
                statistic="mean",
                bins=nshells,
            )[0]
            / stats.binned_statistic(
                pos_prof_n, mass_prof_n, statistic="mean", bins=nshells
            )[0]
        )

    # Calculate bin volumes
    vol_prof_p_new = np.array(
        [
            4 / 3 * np.pi * (bins_p[i + 1] ** 3 - bins_p[i] ** 3)
            for i in range(len(bins_p) - 1)
        ]
    )
    vol_prof_n_new = np.array(
        [
            4 / 3 * np.pi * (bins_n[i + 1] ** 3 - bins_n[i] ** 3)
            for i in range(len(bins_n) - 1)
        ]
    )

    # Sum masses in bins
    mass_prof_p_new = stats.binned_statistic(
        pos_prof_p, mass_prof_p, statistic="sum", bins=nshells
    )[0]
    mass_prof_n_new = stats.binned_statistic(
        pos_prof_n, mass_prof_n, statistic="sum", bins=nshells
    )[0]

    # Calculate densities
    rho_prof_p_new = mass_prof_p_new / vol_prof_p_new
    rho_prof_n_new = mass_prof_n_new / vol_prof_n_new

    # Bin centers
    pos_prof_p_new = np.array(
        [(bins_p[i] + bins_p[i + 1]) / 2 for i in range(len(bins_p) - 1)]
    )
    pos_prof_n_new = np.array(
        [(bins_n[i] + bins_n[i + 1]) / 2 for i in range(len(bins_n) - 1)]
    )

    return (
        pos_prof_p_new,
        pos_prof_n_new,
        vel_prof_p_new,
        vel_prof_n_new,
        rho_prof_p_new,
        rho_prof_n_new,
        mass_prof_p_new,
        mass_prof_n_new,
        xnuc_prof_p_new,
        xnuc_prof_n_new,
    )


def export_profile_to_csvy(
    pos_prof: np.ndarray,
    vel_prof: np.ndarray,
    rho_prof: np.ndarray,
    mass_prof: np.ndarray,
    xnuc_prof: dict,
    time: u.Quantity,
    filename: str | Path,
    nshells: int,
    overwrite: bool = False,
) -> str:
    """
    Export a 1D profile to CSVY format for TARDIS.

    Parameters
    ----------
    pos_prof : numpy.ndarray
        Position profile (cm).
    vel_prof : numpy.ndarray
        Velocity profile (cm/s).
    rho_prof : numpy.ndarray
        Density profile (g/cm^3).
    mass_prof : numpy.ndarray
        Mass profile (g).
    xnuc_prof : dict
        Nuclear fraction profiles keyed by species name.
    time : astropy.units.Quantity
        Time of the snapshot.
    filename : str or pathlib.Path
        Name of the exported file.
    nshells : int
        Number of shells to export.
    overwrite : bool, optional
        If True, will overwrite existing files. Default: False

    Returns
    -------
    str
        Actual filename of the saved file.
    """
    filename = Path(filename)
    species = list(xnuc_prof.keys())

    # Handle filename collision
    if filename.suffix == ".csvy":
        filename = filename.with_suffix("")

    if filename.with_suffix(".csvy").exists() and not overwrite:
        i = 0
        while (
            filename.with_name(f"{filename.stem}_{i}")
            .with_suffix(".csvy")
            .exists()
        ):
            i += 1
        filename = filename.with_name(f"{filename.stem}_{i}")

    filename = filename.with_suffix(".csvy")

    with open(filename, "w") as f:
        # Write YAML header
        f.write(
            "".join(
                [
                    "---\n",
                    "name: csvy_full\n",
                    f"model_density_time_0: {time.to(u.day):g}\n",
                    f"model_isotope_time_0: {time.to(u.day):g}\n",
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

        for spec in species:
            f.write(
                "".join(
                    [
                        f"    -  name: {spec.capitalize()}\n",
                        f"       desc: fractional {spec.capitalize()} abundance.\n",
                    ]
                )
            )

        f.write("\n---\n")

        # Write CSV header
        datastring = ["velocity,", "density,"]
        for spec in species[:-1]:
            datastring.append(f"{spec.capitalize()},")
        datastring.append(f"{species[-1].capitalize()}")
        f.write("".join(datastring))

        # Prepare data arrays
        exp = [vel_prof, rho_prof]
        for spec in species:
            exp.append(xnuc_prof[spec])

        # Write data rows
        inds = np.linspace(0, len(exp[0]) - 1, num=nshells, dtype=int)
        for i in inds:
            f.write("\n")
            for ii in range(len(exp) - 1):
                f.write(f"{exp[ii][i]:g},")
            f.write(f"{exp[-1][i]:g}")

    return str(filename)


def plot_profile(
    pos_prof_p: np.ndarray,
    pos_prof_n: np.ndarray,
    vel_prof_p: np.ndarray,
    vel_prof_n: np.ndarray,
    rho_prof_p: np.ndarray,
    rho_prof_n: np.ndarray,
    xnuc_prof_p: dict,
    xnuc_prof_n: dict,
    time: u.Quantity,
    save: str | Path | None = None,
    dpi: int = 600,
    **kwargs,
) -> "Figure":
    """
    Plot 1D profiles for both positive and negative directions.

    Parameters
    ----------
    pos_prof_p : numpy.ndarray
        Position profile in positive direction.
    pos_prof_n : numpy.ndarray
        Position profile in negative direction.
    vel_prof_p : numpy.ndarray
        Velocity profile in positive direction.
    vel_prof_n : numpy.ndarray
        Velocity profile in negative direction.
    rho_prof_p : numpy.ndarray
        Density profile in positive direction.
    rho_prof_n : numpy.ndarray
        Density profile in negative direction.
    xnuc_prof_p : dict
        Nuclear fraction profiles in positive direction.
    xnuc_prof_n : dict
        Nuclear fraction profiles in negative direction.
    time : astropy.units.Quantity
        Time of the snapshot.
    save : str or pathlib.Path, optional
        Path to save the figure. Default: None
    dpi : int, optional
        DPI for saved figure. Default: 600
    **kwargs
        Additional keyword arguments passed to matplotlib.pyplot.plot()

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure object.
    """
    species = list(xnuc_prof_p.keys())
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=[9.8, 9.6])

    # Positive direction plots
    ax1.plot(
        pos_prof_p,
        rho_prof_p / max(rho_prof_p),
        label="Density",
        **kwargs,
    )
    ax1.plot(
        pos_prof_p,
        vel_prof_p / max(vel_prof_p),
        label="Velocity",
        **kwargs,
    )
    for spec in species:
        ax1.plot(
            pos_prof_p,
            xnuc_prof_p[spec],
            label=spec.capitalize(),
            **kwargs,
        )

    ax1.grid()
    ax1.set_ylabel("Profile (arb. unit)")
    ax1.set_title("Profiles along the positive axis")

    # Negative direction plots
    ax2.plot(
        pos_prof_n,
        rho_prof_n / max(rho_prof_n),
        label="Density",
        **kwargs,
    )
    ax2.plot(
        pos_prof_n,
        vel_prof_n / max(vel_prof_n),
        label="Velocity",
        **kwargs,
    )
    for spec in species:
        ax2.plot(
            pos_prof_n,
            xnuc_prof_n[spec],
            label=spec.capitalize(),
            **kwargs,
        )

    ax2.grid()
    ax2.set_ylabel("Profile (arb. unit)")
    ax2.set_xlabel("Radial position (cm)")
    ax2.set_title("Profiles along the negative axis")

    # Styling
    fig.tight_layout()

    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(1.05, 1.05),
        title=f"Time = {time}",
    )

    if save is not None:
        plt.savefig(save, bbox_inches="tight", dpi=dpi)

    return fig
