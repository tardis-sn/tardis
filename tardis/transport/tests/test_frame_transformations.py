import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_allclose

import tardis.transport.frame_transformations as frame_transformations
from tardis import constants as const
from tardis.transport.montecarlo.packets.radiative_packet import RPacket

C_SPEED_OF_LIGHT = const.c.to("cm/s").value
# 13 days [s]; any positive value works because only beta = r / (c t) enters.
TIME_EXPLOSION = 13.0 * u.day.to("s")
# Direction cosines and frequency ratios are O(1) and each transformation is a
# handful of floating-point operations, so round trips are exact to a few ulp
# (~1e-16). 1e-14 allows that rounding while rejecting any algebraic error.
ROUND_TRIP_TOL = 1e-14


def make_packet_at_beta(beta: float, mu: float = 0.0) -> RPacket:
    """
    Build a packet at the radius whose homologous velocity is beta * c.

    The aberration functions read beta = r / (c t) from the packet, so the
    radius is chosen to give the requested beta at TIME_EXPLOSION.
    """
    r = beta * C_SPEED_OF_LIGHT * TIME_EXPLOSION
    return RPacket(r=r, mu=mu, nu=6.0e14, energy=1.0, seed=1963)


# beta = 0.03 (9000 km/s) is typical of supernova ejecta; 0.3 stresses the
# relativistic terms.
@pytest.mark.parametrize("beta", [0.03, 0.3])
@pytest.mark.parametrize("mu", [-1.0, -0.3, 0.0, 0.3, 1.0])
def test_angle_aberration_round_trip(mu: float, beta: float) -> None:
    """
    Claim: lab-to-comoving and comoving-to-lab angle aberration are inverses.

    Regime: special-relativistic aberration at the packet's homologous
    velocity.
    Verification: composition of the two transforms must be the identity on
    [-1, 1] (metamorphic relation), in both orders.
    """
    packet = make_packet_at_beta(beta)

    mu_cmf = frame_transformations.angle_aberration_LF_to_CMF(
        packet, TIME_EXPLOSION, mu
    )
    mu_lab = frame_transformations.angle_aberration_CMF_to_LF(
        packet, TIME_EXPLOSION, mu
    )

    assert_allclose(
        frame_transformations.angle_aberration_CMF_to_LF(
            packet, TIME_EXPLOSION, mu_cmf
        ),
        mu,
        atol=ROUND_TRIP_TOL,
    )
    assert_allclose(
        frame_transformations.angle_aberration_LF_to_CMF(
            packet, TIME_EXPLOSION, mu_lab
        ),
        mu,
        atol=ROUND_TRIP_TOL,
    )


@pytest.mark.parametrize(
    ["beta", "mu_cmf", "expected_mu_lab"],
    [
        # Radial directions are fixed points of aberration.
        (0.3, 1.0, 1.0),
        (0.3, -1.0, -1.0),
        # Emission perpendicular to the flow in the comoving frame is beamed
        # forward to mu_lab = beta.
        (0.3, 0.0, 0.3),
        # At rest (beta = 0, i.e. r = 0) the frames coincide.
        (0.0, 0.42, 0.42),
    ],
)
def test_angle_aberration_cmf_to_lf_limits(
    beta: float, mu_cmf: float, expected_mu_lab: float
) -> None:
    """
    Claim: comoving-to-lab aberration reproduces the analytic limits of
    mu_lab = (mu_cmf + beta) / (1 + beta * mu_cmf).

    Regime: radial directions, perpendicular comoving emission, and the
    zero-velocity limit.
    Verification: closed-form values of the aberration formula.
    """
    packet = make_packet_at_beta(beta)

    mu_lab = frame_transformations.angle_aberration_CMF_to_LF(
        packet, TIME_EXPLOSION, mu_cmf
    )

    assert_allclose(mu_lab, expected_mu_lab, atol=ROUND_TRIP_TOL)


@pytest.mark.parametrize("beta", [0.03, 0.3])
@pytest.mark.parametrize("mu_lab", [-1.0, -0.3, 0.0, 0.3, 1.0])
def test_doppler_factor_partial_relativity_inverse(
    mu_lab: float, beta: float
) -> None:
    """
    Claim: in partial relativity the inverse Doppler factor is the reciprocal
    of the Doppler factor at the same lab-frame angle.

    Regime: first-order (partial) relativity, D = 1 - beta * mu_lab.
    Verification: D * D_inverse == 1, so a lab -> comoving -> lab frequency
    transformation returns the original frequency.
    """
    velocity = beta * C_SPEED_OF_LIGHT

    doppler_factor = frame_transformations.get_doppler_factor(
        velocity, mu_lab, False
    )
    inverse_doppler_factor = frame_transformations.get_inverse_doppler_factor(
        velocity, mu_lab, False
    )

    assert_allclose(
        doppler_factor * inverse_doppler_factor, 1.0, rtol=ROUND_TRIP_TOL
    )


@pytest.mark.parametrize("beta", [0.03, 0.3])
@pytest.mark.parametrize("mu_lab", [-1.0, -0.3, 0.0, 0.3, 1.0])
def test_doppler_factor_full_relativity_inverse(
    mu_lab: float, beta: float
) -> None:
    """
    Claim: in full relativity a lab -> comoving -> lab frequency round trip
    returns the original frequency.

    Regime: full special relativity. The forward factor
    gamma * (1 - beta * mu_lab) takes the lab-frame angle, while the inverse
    gamma * (1 + beta * mu_cmf) takes the comoving-frame angle, so the
    direction must be aberrated between the two.
    Verification: gamma^2 (1 - beta mu)(1 + beta mu') = 1 when
    mu' = (mu - beta) / (1 - beta mu), the exact relativistic aberration.
    """
    velocity = beta * C_SPEED_OF_LIGHT
    packet = make_packet_at_beta(beta)
    mu_cmf = frame_transformations.angle_aberration_LF_to_CMF(
        packet, TIME_EXPLOSION, mu_lab
    )

    doppler_factor = frame_transformations.get_doppler_factor(
        velocity, mu_lab, True
    )
    inverse_doppler_factor = frame_transformations.get_inverse_doppler_factor(
        velocity, mu_cmf, True
    )

    assert_allclose(
        doppler_factor * inverse_doppler_factor, 1.0, rtol=ROUND_TRIP_TOL
    )


@pytest.mark.parametrize(
    "distance_fraction",
    # Distance travelled as a fraction of c * t: zero (energy at the current
    # position), and 1e-3, a typical line-resonance distance at a 0.1%
    # comoving frequency offset.
    [0.0, 1e-3],
)
@pytest.mark.parametrize("mu", [-1.0, -0.3, 0.0, 0.3, 1.0])
def test_calc_packet_energy_is_comoving_energy_after_move(
    mu: float, distance_fraction: float
) -> None:
    """
    Claim: calc_packet_energy returns the packet's comoving-frame energy at the
    point reached after travelling ``distance_trace``, to first order in beta.

    Regime: partial relativity, homologous flow, beta = 0.03.
    Verification: the new position follows from the law of cosines, and the
    energy is E_lab * (1 - beta_new * mu_new) with beta_new = r_new / (c t),
    the first-order Doppler factor evaluated there.
    """
    beta = 0.03
    energy_lab = 0.9
    ct = C_SPEED_OF_LIGHT * TIME_EXPLOSION
    distance = distance_fraction * ct
    packet = make_packet_at_beta(beta, mu=mu)
    packet.energy = energy_lab
    r = packet.r

    r_new = np.sqrt(r * r + distance * distance + 2.0 * r * distance * mu)
    mu_new = (r * mu + distance) / r_new
    expected_energy = energy_lab * (1.0 - (r_new / ct) * mu_new)

    energy_cmf = frame_transformations.calc_packet_energy(
        packet, distance, TIME_EXPLOSION
    )

    assert_allclose(energy_cmf, expected_energy, rtol=ROUND_TRIP_TOL)
