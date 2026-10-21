from collections.abc import Callable

import numpy as np
import numpy.testing as npt
import pytest

from tardis.energy_input.transport.gamma_ray_interactions import (
    compton_scatter,
    get_compton_angle,
    pair_creation_packet,
    scatter_type,
)
from tardis.energy_input.transport.GXPacket import GXPacket, GXPacketStatus
from tardis.energy_input.util import (
    ELECTRON_MASS_ENERGY_KEV,
    H_CGS_KEV,
    angle_aberration_gamma,
)


def test_get_compton_angle() -> None:
    """Verify that sampled energies remain within the Compton edge.

    Claim: A photon cannot gain energy from a stationary electron or lose more
    than the energy transferred in a 180-degree backscatter.
    Regime: 1000 keV photons scattering from stationary electrons.
    Verification: Analytic forward- and backscatter energy limits.
    """
    initial_energy = 1000.0  # keV
    minimum_scattered_energy = initial_energy / (
        1.0 + 2.0 * initial_energy / ELECTRON_MASS_ENERGY_KEV
    )
    roundoff_tolerance = 10.0 * np.finfo(np.float64).eps * initial_energy

    for _ in range(100):
        angle, lost_energy, scattered_energy = get_compton_angle(initial_energy)

        assert 0.0 <= angle <= np.pi
        assert minimum_scattered_energy - roundoff_tolerance <= scattered_energy
        assert scattered_energy <= initial_energy + roundoff_tolerance
        assert -roundoff_tolerance <= lost_energy
        assert lost_energy <= (
            initial_energy - minimum_scattered_energy + roundoff_tolerance
        )


def test_get_compton_angle_becomes_forward_peaked(
    set_seed_fixture: Callable[[int], None],
) -> None:
    """Verify stronger forward scattering at higher photon energy.

    Claim: Increasing photon energy makes Klein-Nishina scattering more
    forward-peaked.
    Regime: Compare 1 keV and 1000 keV photons using paired random quantiles.
    Verification: Metamorphic ordering of the mean scattering cosine.
    """
    sample_size = 4096
    seed = 1963

    set_seed_fixture(seed)
    low_energy_cosines = np.array(
        [np.cos(get_compton_angle(1.0)[0]) for _ in range(sample_size)]
    )
    set_seed_fixture(seed)
    high_energy_cosines = np.array(
        [np.cos(get_compton_angle(1000.0)[0]) for _ in range(sample_size)]
    )
    paired_difference = high_energy_cosines - low_energy_cosines
    standard_error = np.std(paired_difference, ddof=1) / np.sqrt(sample_size)

    assert np.mean(paired_difference) > 5.0 * standard_error, (
        "Klein-Nishina sampling did not become significantly more "
        "forward-peaked with photon energy"
    )


def test_compton_scatter(basic_gamma_ray: GXPacket) -> None:
    """Verify that scattering applies the requested polar angle.

    Claim: The angle between the incident and scattered directions equals the
    requested Compton angle.
    Regime: A sub-relativistic homologous flow with v/c approximately 0.005.
    Verification: Transform the returned lab-frame direction independently back
    to the comoving frame and apply the geometric dot-product identity.
    """
    basic_gamma_ray.location = np.array([3.0e13, 0.0, 0.0])
    basic_gamma_ray.direction = np.array([0.8, 0.6, 0.0])
    basic_gamma_ray.time_start = 2.0e5
    compton_angle = 0.5  # radians
    initial_comoving_direction = angle_aberration_gamma(
        basic_gamma_ray.direction,
        basic_gamma_ray.location,
        basic_gamma_ray.time_start,
    )

    scattered_direction = compton_scatter(basic_gamma_ray, compton_angle)
    scattered_comoving_direction = angle_aberration_gamma(
        scattered_direction,
        basic_gamma_ray.location,
        basic_gamma_ray.time_start,
    )

    npt.assert_allclose(
        np.dot(initial_comoving_direction, scattered_comoving_direction),
        np.cos(compton_angle),
        rtol=1.0e-12,
        err_msg="Scattered direction does not have the requested polar angle",
    )


def test_compton_scatter_preserves_unit_direction(
    basic_gamma_ray: GXPacket,
) -> None:
    """Verify that Compton scattering returns a unit direction.

    Claim: Rotating a photon direction does not change its vector norm.
    Regime: A sub-relativistic homologous flow with v/c approximately 0.005.
    Verification: Direction-vector norm invariance after transforming the
    scattered direction back to the comoving frame.
    """
    basic_gamma_ray.location = np.array([3.0e13, 0.0, 0.0])
    basic_gamma_ray.direction = np.array([0.8, 0.6, 0.0])
    basic_gamma_ray.time_start = 2.0e5

    scattered_direction = compton_scatter(basic_gamma_ray, 0.5)
    scattered_comoving_direction = angle_aberration_gamma(
        scattered_direction,
        basic_gamma_ray.location,
        basic_gamma_ray.time_start,
    )

    npt.assert_allclose(
        np.linalg.norm(scattered_comoving_direction),
        1.0,
        rtol=1.0e-12,
        err_msg="Compton scattering changed the direction-vector norm",
    )


def test_pair_creation(basic_gamma_ray):
    """

    Parameters
    ----------
    basic_gamma_ray : GammaRay object
    """
    np.random.seed(2)

    initial_direction = basic_gamma_ray.direction
    basic_gamma_ray.nu_cmf = 2 * ELECTRON_MASS_ENERGY_KEV / H_CGS_KEV

    pair_creation_packet(basic_gamma_ray)

    npt.assert_almost_equal(
        basic_gamma_ray.nu_cmf, ELECTRON_MASS_ENERGY_KEV / H_CGS_KEV
    )
    for i, j in zip(initial_direction, basic_gamma_ray.direction, strict=True):
        assert i != j


@pytest.mark.parametrize(
    ["compton_opacity", "photoabsorption_opacity", "total_opacity", "expected"],
    [
        (1, 0, 1, GXPacketStatus.COMPTON_SCATTER),
        (0, 1, 1, GXPacketStatus.PHOTOABSORPTION),
        (0, 0, 1, GXPacketStatus.PAIR_CREATION),
    ],
)
def test_scatter_type(
    compton_opacity, photoabsorption_opacity, total_opacity, expected
):
    """Test the scattering type

    Parameters
    ----------
    compton_opacity : float
    photoabsorption_opacity : float
    total_opacity : float
    expected : list
        Expected parameters
    """
    actual = scatter_type(
        compton_opacity, photoabsorption_opacity, total_opacity
    )
    assert actual == expected
