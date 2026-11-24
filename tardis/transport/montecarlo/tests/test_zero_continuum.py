"""Eventful zero-continuum limits of the current classic and IIP kernels."""

import numpy as np
import numpy.testing as npt
import pytest
from numpy.typing import NDArray
from scipy.stats import t

from tardis.conftest import set_numba_seed
from tardis.model.geometry.radial1d_homologous import (
    NumbaHomologousRadial1DGeometry,
)
from tardis.opacities.continuum.continuum_state_numba import (
    ContinuumOpacityStateNumba,
)
from tardis.opacities.opacities import FF_OPAC_CONST, K_B, H
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
)
from tardis.transport.montecarlo.configuration.constants import SIGMA_THOMSON
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    init_estimators_bulk,
)
from tardis.transport.montecarlo.estimators.estimators_continuum import (
    init_estimators_continuum,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    init_estimators_line,
)
from tardis.transport.montecarlo.interaction_events import LineInteractionType
from tardis.transport.montecarlo.macro_atom import MacroAtomTransitionType
from tardis.transport.montecarlo.modes.classic.packet_propagation import (
    packet_propagation as propagate_classic,
)
from tardis.transport.montecarlo.modes.iip.packet_propagation import (
    packet_propagation as propagate_iip,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    PacketStatus,
    RPacket,
)
from tardis.transport.montecarlo.packets.trackers.tracker_full import (
    TrackerFull,
)

# Predeclared design, shared by the zero-opacity comparison and positive control.
BATCH_COUNT = 32
PACKETS_PER_BATCH = 128
OBSERVABLE_NAMES = (
    "escaped energy",
    "escaped energy-frequency moment",
    "bulk path energy",
    "bulk path energy-frequency moment",
    "blue-line intensity",
)
FAMILY_ALPHA = 0.001
SIGNIFICANCE_THRESHOLD = t.isf(
    FAMILY_ALPHA / (2 * len(OBSERVABLE_NAMES)), BATCH_COUNT - 1
)


@pytest.fixture(scope="module")
def zero_continuum_transport_batches() -> NDArray[np.float64]:
    """Sample matched radiative macro atoms with independent packet seed batches."""
    # cm, s, Hz, erg, K. Thomson depth across the ejecta is 0.9.
    time_explosion = 1.0e6
    r_inner = np.array([7.0e14, 8.0e14])
    r_outer = np.array([8.0e14, 1.0e15])
    geometry = NumbaHomologousRadial1DGeometry(
        r_inner,
        r_outer,
        r_inner / time_explosion,
        r_outer / time_explosion,
        time_explosion,
    )
    source_nu = 4.0e14
    source_energy = 1.0
    temperature = 1.0e4
    opacity = OpacityStateNumba(
        np.full(2, 3.0e-15 / SIGMA_THOMSON),
        np.full(2, temperature),
        # Reachable absorption line, lower radiative deactivation/sentinel line.
        np.array([3.98e14, 2.0e14]),
        np.array([[2.0, 2.0], [0.0, 0.0]]),
        np.ones((1, 2)),
        np.zeros(2, dtype=np.int64),
        np.array([0, 1], dtype=np.int64),
        np.array([MacroAtomTransitionType.BB_EMISSION], dtype=np.int64),
        np.zeros(1, dtype=np.int64),
        np.ones(1, dtype=np.int64),
    )
    # No active BF continua: an active zero-filled interpolation grid would
    # divide by zero. Positive electron opacity keeps chi_e/chi_total valid.
    continuum = ContinuumOpacityStateNumba(
        np.empty(0),
        np.empty((0, 2)),
        np.empty(0),
        np.empty(0),
        np.array([0], dtype=np.int64),
        np.empty((0, 2)),
        np.empty(0),
        np.empty(0),
        np.zeros(2),
        np.empty((0, 2)),
        np.empty(0, dtype=np.int64),
        0,
        # IIP's absorbing-chain jump is exactly to the same radiative state;
        # classic samples the sole BB channel directly. Draw counts differ.
        np.ones((2, 1, 1)),
        # One macro level, so global and reduced level indices coincide.
        np.arange(1, dtype=np.int64),
        np.arange(1, dtype=np.int64),
    )
    configuration = MonteCarloConfiguration()
    configuration.ENABLE_FULL_RELATIVITY = True
    configuration.LINE_INTERACTION_TYPE = LineInteractionType.MACROATOM

    # Axis 0: classic off, IIP zero continuum, IIP finite free-free control.
    observables = np.zeros((3, BATCH_COUNT, len(OBSERVABLE_NAMES)))
    for mode_idx in range(3):
        # The control has chi_ff=3e-15 cm^-1 at the source frequency, with
        # real FF absorption followed by the same forced radiative BB channel.
        continuum.ff_opacity_factor[:] = (
            3.0e-15
            * source_nu**3
            / (
                FF_OPAC_CONST
                * (1 - np.exp(-H * source_nu / (K_B * temperature)))
            )
            if mode_idx == 2
            else 0.0
        )
        for batch_idx in range(BATCH_COUNT):
            bulk = init_estimators_bulk(2)
            line = init_estimators_line((2, 2))
            continuous = init_estimators_continuum((0, 2), 2)
            for packet_idx in range(PACKETS_PER_BATCH):
                # Disjoint seeds define independent batches; resetting the
                # same packet seed across modes retains any useful pairing.
                seed = 104729 + batch_idx * PACKETS_PER_BATCH + packet_idx
                packet = RPacket(
                    r=7.5e14,
                    mu=0.6,
                    nu=source_nu,
                    energy=source_energy,
                    seed=seed,
                    index=packet_idx,
                )
                tracker = TrackerFull(32, 2)
                set_numba_seed(seed)
                if mode_idx == 0:
                    propagate_classic(
                        packet,
                        geometry,
                        time_explosion,
                        opacity,
                        bulk,
                        line,
                        tracker,
                        configuration,
                    )
                else:
                    propagate_iip(
                        packet,
                        geometry,
                        time_explosion,
                        opacity,
                        continuum,
                        bulk,
                        line,
                        continuous,
                        tracker,
                        configuration,
                    )
                if packet.status == PacketStatus.EMITTED:
                    observables[mode_idx, batch_idx, 0] += packet.energy
                    observables[mode_idx, batch_idx, 1] += (
                        packet.energy * packet.nu / source_nu
                    )
            # Fixed units/scales merely nondimensionalize the batch means.
            observables[mode_idx, batch_idx, 2:] = (
                bulk.mean_intensity_total.sum() / 1.0e14,
                bulk.mean_frequency.sum() / (1.0e14 * source_nu),
                line.mean_intensity_blueward.sum() * source_nu,
            )
            observables[mode_idx, batch_idx] /= (
                PACKETS_PER_BATCH * source_energy
            )
    return observables


def test_zero_continuum_reduces_to_line_only_transport(
    zero_continuum_transport_batches: NDArray[np.float64],
) -> None:
    """Verify the eventful zero-continuum limiting distribution (verification).

    Claim: zero BF/FF opacity reduces FULL/MACROATOM IIP transport to classic
    line/electron transport in distribution. Regime: homologous finite-velocity
    ejecta, positive Thomson opacity, one reachable tau=2 line with exclusively
    radiative deactivation into the same lower line in both representations.
    Verification: the independent classic propagator supplies the limiting-case
    oracle; no packet states or event histories are compared across modes.

    Before sampling, fix B=32 independent batches of 128 packets and the five
    observables above. For observable j and batch b, D_bj=IIP_bj-classic_bj.
    SE_j=sqrt(sum_b(D_bj-mean(D_j))**2/[B*(B-1)]). Equivalently,
    Var(D)=Var(IIP)+Var(classic)-2*Cov(IIP,classic): shared seeds induce pairing
    even though split/Markov draws destroy trajectory equality. Using the
    measured differences includes that correlation without assuming it zero.
    Each batch mean has 128 independent histories; a Student-t bound with
    B-1 degrees of freedom approximates its sampling distribution. Bonferroni
    alpha/(2*5) controls the two-sided family error at 0.001 regardless of
    correlations among observables. The fixed bound is t_31(0.9999), about 4.22
    standard errors, with no relative-tolerance tuning or absolute noise floor.
    At exactly zero SE, exact mean agreement has standardized difference zero;
    a nonzero mean has signed infinity and fails the same significance bound.

    A finite-FF actual-path control must exceed the SAME zero-limit MC bound
    in the preselected escaped energy-frequency moment. This demonstrates
    resolution of continuum contamination/process handoff, not just agreement
    of noisy samples. Detects residual continuum opacity or lost line/electron
    handling; equivalent deactivation channels isolate transport from atom data.
    """
    observables = zero_continuum_transport_batches
    differences = observables[1] - observables[0]
    mean_differences = differences.mean(axis=0)
    standard_errors = differences.std(axis=0, ddof=1) / np.sqrt(BATCH_COUNT)
    standardized_differences = np.copysign(
        np.full_like(mean_differences, np.inf), mean_differences
    )
    np.divide(
        mean_differences,
        standard_errors,
        out=standardized_differences,
        where=standard_errors != 0.0,
    )
    standardized_differences[
        (standard_errors == 0.0) & (mean_differences == 0.0)
    ] = 0.0
    npt.assert_array_less(
        np.abs(standardized_differences),
        np.full(len(OBSERVABLE_NAMES), SIGNIFICANCE_THRESHOLD),
        err_msg=f"Zero-continuum reduction failed for {OBSERVABLE_NAMES}",
    )
    control_difference = (observables[2] - observables[0]).mean(axis=0)
    assert (
        abs(control_difference[1]) > SIGNIFICANCE_THRESHOLD * standard_errors[1]
    ), (
        "Finite free-free opacity must shift the escaped frequency moment "
        "beyond the same predeclared zero-continuum noise threshold"
    )
