from pathlib import Path

import numpy as np
from astropy import units as u

from tardis.io.model.arepo.data import ArepoData


def read_arepo_snapshot(
    filename: str | Path,
    species: list[str],
    speciesfile: str | Path,
    alpha: float = 0.0,
    beta: float = 0.0,
    gamma: float = 0.0,
) -> ArepoData:
    """
    Read an Arepo snapshot file and return ArepoData.

    Loads relevant data from an Arepo snapshot for conversion to TARDIS model.
    Requires arepo-snap-util to be installed. The snapshot is rotated according
    to the provided Euler angles.

    Parameters
    ----------
    filename : str or pathlib.Path
        Path to the Arepo snapshot file.
    species : list of str
        Names of the species to be exported. Must match the species-file
        of the Arepo simulation.
    speciesfile : str or pathlib.Path
        File specifying the species used in the Arepo simulation.
    alpha : float, optional
        Euler angle alpha (yaw) for rotation of the desired line-of-sight
        to the x-axis. Default: 0.0
    beta : float, optional
        Euler angle beta (pitch) for rotation of the desired line-of-sight
        to the x-axis. Default: 0.0
    gamma : float, optional
        Euler angle gamma (roll) for rotation of the desired line-of-sight
        to the x-axis. Default: 0.0

    Returns
    -------
    ArepoData
        Data structure containing the Arepo snapshot data.

    Raises
    ------
    ImportError
        If arepo-snap-util is not installed.
    """
    try:
        import gadget_snap
    except ModuleNotFoundError:
        raise ImportError(
            "Please make sure you have arepo-snap-util installed if you want to directly import Arepo snapshots."
        )

    filename = Path(filename)
    speciesfile = Path(speciesfile)

    # Read species file and find indices
    species_full = np.genfromtxt(speciesfile, skip_header=1, dtype=str).T[0]
    spec_ind = []
    for spec in species:
        spec_ind.append(np.where(species_full == spec)[0][0])
    spec_ind = np.array(spec_ind)

    # Load snapshot
    s = gadget_snap.gadget_snapshot(
        str(filename),
        hdf5=True,
        quiet=True,
        lazy_load=True,
        loadonlytype=[0],
    )

    # Compute rotation matrix from Euler angles
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

    s.rotateto(rotmat[0], dir2=rotmat[1], dir3=rotmat[2])

    # Extract data
    time = u.Quantity(s.time, "s")
    pos = np.array(s.data["pos"]).T
    # Update position to CoM frame
    for i in range(3):
        pos[i] -= s.centerofmass()[i]
    rho = np.array(s.data["rho"])
    mass = np.array(s.data["mass"])
    vel = np.array(s.data["vel"]).T

    # Extract nuclear fractions
    nuc_dict = {}
    for i, spec in enumerate(species):
        nuc_dict[spec] = np.array(s.data["xnuc"][:, spec_ind[i]])

    return ArepoData(
        time=time,
        position=pos,
        velocities=vel,
        densities=rho,
        mass=mass,
        isotope_dict=nuc_dict,
    )
