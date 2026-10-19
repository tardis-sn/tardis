"""Eventful zero-continuum limits of the current classic and IIP kernels."""

from collections.abc import Callable

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
    EstimatorsLine,
    init_estimators_line,
)
from tardis.transport.montecarlo.interaction_events import LineInteractionType
from tardis.transport.montecarlo.macro_atom import MacroAtomTransitionType
from tardis.transport.montecarlo.modes.classic.packet_propagation import (
    packet_propagation as propagate_classic,
)
from tardis.transport.montecarlo.modes.homologous_rad_packet_transport import (
    trace_packet,
)
from tardis.transport.montecarlo.modes.iip.packet_propagation import (
    packet_propagation as propagate_iip,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    PacketStatus,
    RPacket,
)
from tardis.transport.montecarlo.packets.trackers.tracker_full import (
    TrackerFull,
)
from tardis.transport.montecarlo.utils import get_random_mu

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


def calculate_standardized_differences(
    mean_differences: NDArray[np.float64],
    standard_errors: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Keep exact agreement at zero SE and reject any resolved constant bias."""
    standardized = np.copysign(
        np.full_like(mean_differences, np.inf), mean_differences
    )
    np.divide(
        mean_differences,
        standard_errors,
        out=standardized,
        where=standard_errors != 0.0,
    )
    standardized[(standard_errors == 0.0) & (mean_differences == 0.0)] = 0.0
    return standardized


@pytest.mark.parametrize(
    "batch_differences,expected",
    [
        ([0.0, 0.0], 0.0),
        ([1.0, 1.0], np.inf),
        ([-1.0, -1.0], -np.inf),
        ([1.0, 3.0], 2.0),
    ],
)
def test_standardized_differences_include_exact_paired_agreement(
    batch_differences: list[float],
    expected: float,
) -> None:
    """Verify exact agreement/bias at zero noise and the ordinary finite ratio."""
    differences = np.array(batch_differences)
    actual = calculate_standardized_differences(
        np.array([differences.mean()]),
        np.array([differences.std(ddof=1) / np.sqrt(len(differences))]),
    )
    npt.assert_array_equal(actual, [expected])


@pytest.fixture(scope="module")
def zero_continuum_transport_batches() -> dict[str, NDArray]:
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
    event_counts = np.zeros((3, 3), dtype=np.int64)
    ff_heating = np.zeros(3)
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
                for event_idx, event in enumerate(
                    (
                        InteractionType.LINE,
                        InteractionType.ESCATTERING,
                        InteractionType.CONTINUUM_PROCESS,
                    )
                ):
                    event_counts[mode_idx, event_idx] += np.count_nonzero(
                        tracker.interaction_type[: tracker.event_id] == event
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
            ff_heating[mode_idx] += continuous.ff_heating_estimator.sum()
    return {
        "observables": observables,
        "event_counts": event_counts,
        "ff_heating": ff_heating,
    }


def test_zero_continuum_reduces_to_line_only_transport(
    zero_continuum_transport_batches: dict[str, NDArray],
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
    observables = zero_continuum_transport_batches["observables"]
    differences = observables[1] - observables[0]
    standard_errors = differences.std(axis=0, ddof=1) / np.sqrt(BATCH_COUNT)
    standardized_differences = calculate_standardized_differences(
        differences.mean(axis=0), standard_errors
    )
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


def test_zero_continuum_batches_exercise_current_event_paths(
    zero_continuum_transport_batches: dict[str, NDArray],
) -> None:
    """Audit sampled process coverage separately from physical verification."""
    counts = zero_continuum_transport_batches["event_counts"]
    assert np.all(counts[:, :2] > 0), (
        "Every mode must exercise lines and electrons"
    )
    npt.assert_array_equal(counts[:2, 2], 0)
    assert counts[2, 2] > 0, (
        "The finite-FF control must traverse continuum_event"
    )
    heating = zero_continuum_transport_batches["ff_heating"]
    npt.assert_array_equal(heating[:2], 0.0)
    assert heating[2] > 0.0


def test_zero_continuum_trace_consumes_extra_split_draw(
    parametrized_packet: RPacket,
    homologous_radial_1d_geometry: NumbaHomologousRadial1DGeometry,
    classic_opacity_state: OpacityStateNumba,
    line_estimators: EstimatorsLine,
    set_seed_fixture: Callable[[int], None],
) -> None:
    """Characterize identical next events with different subsequent RNG draws.

    Regime: FULL homologous trace, positive continuous opacity, escat_prob=1,
    electron distance before every line/boundary. Reset seed for EACH call.
    Verification: analytic exponential distance and an independent RandomState
    sequence show one optical-depth draw off versus depth plus split draw on.
    Geometry tolerance rtol=1e-12 follows the existing trace tests; the RNG
    conversion allows roundoff only. Removing the split or changing the chosen
    interaction falsifies this diagnostic. This is not RNG equivalence or an
    assertion that transported packet histories agree.
    """
    seed = 1963
    draws = np.random.RandomState(seed).random_sample(3)
    opacity = 1.0e-10  # cm^-1: mean free path 1e10 cm, boundary >1e13 cm
    # Safe retained redward lines: no resonance before the electron event.
    classic_opacity_state.line_list_nu[:] = [3.0e14, 2.0e14]
    results = []
    following_draws = []
    for enabled in (False, True):
        parametrized_packet.next_line_id = 0
        set_seed_fixture(seed)
        results.append(
            trace_packet(
                parametrized_packet,
                homologous_radial_1d_geometry,
                homologous_radial_1d_geometry.time_explosion,
                classic_opacity_state,
                line_estimators,
                opacity,
                1.0,
                enabled,
                True,
                False,
            )
        )
        # Read the very next draw in the same Numba RNG through existing code.
        following_draws.append((get_random_mu() + 1.0) / 2.0)
    npt.assert_allclose(
        [result[0] for result in results],
        -np.log(draws[0]) / opacity,
        rtol=1.0e-12,
        atol=0.0,
    )
    npt.assert_array_equal(
        [result[1:] for result in results],
        [(InteractionType.ESCATTERING, 1)] * 2,
    )
    npt.assert_allclose(following_draws, draws[1:], rtol=0.0, atol=2.0e-16)
    assert following_draws[0] != following_draws[1]
