from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numba import boolean, int64
from numba.experimental import jitclass

from tardis.transport.montecarlo.configuration import montecarlo_globals
from tardis.transport.montecarlo.interaction_events import (
    LineInteractionType,
)

if TYPE_CHECKING:
    from tardis.transport.montecarlo.modes.classic.solver import (
        MCTransportSolverClassic,
    )
    from tardis.transport.montecarlo.modes.iip.solver import (
        MCTransportSolverIIP,
    )
    from tardis.transport.montecarlo.modes.nonhomologous.solver import (
        MCTransportSolverNonhomologous,
    )

numba_config_spec = [
    ("ENABLE_FULL_RELATIVITY", boolean),
    ("MONTECARLO_SEED", int64),
    ("LINE_INTERACTION_TYPE", int64),
    ("PACKET_SEEDS", int64[:]),
    ("DISABLE_ELECTRON_SCATTERING", boolean),
    ("DISABLE_LINE_SCATTERING", boolean),
    ("INITIAL_TRACKING_ARRAY_LENGTH", int64),
    ("LEGACY_MODE_ENABLED", boolean),
]


@jitclass(numba_config_spec)
class MonteCarloConfiguration:
    def __init__(self):
        self.ENABLE_FULL_RELATIVITY = False
        self.MONTECARLO_SEED = 0
        self.LINE_INTERACTION_TYPE = 0
        self.PACKET_SEEDS = np.empty(1, dtype=np.int64)
        self.DISABLE_ELECTRON_SCATTERING = False
        self.DISABLE_LINE_SCATTERING = False

        self.INITIAL_TRACKING_ARRAY_LENGTH = 0
        self.LEGACY_MODE_ENABLED = False


def configuration_initialize(
    config: MonteCarloConfiguration,
    transport: (
        MCTransportSolverClassic
        | MCTransportSolverNonhomologous
        | MCTransportSolverIIP
    ),
) -> None:
    if transport.line_interaction_type == "macroatom":
        config.LINE_INTERACTION_TYPE = LineInteractionType.MACROATOM
    elif transport.line_interaction_type == "downbranch":
        config.LINE_INTERACTION_TYPE = LineInteractionType.DOWNBRANCH
    elif transport.line_interaction_type == "scatter":
        config.LINE_INTERACTION_TYPE = LineInteractionType.SCATTER
    else:
        raise ValueError(
            f'Line interaction type must be one of "macroatom",'
            f'"downbranch", or "scatter" but is '
            f"{transport.line_interaction_type}"
        )
    config.ENABLE_FULL_RELATIVITY = transport.enable_full_relativity
    config.MONTECARLO_SEED = transport.packet_source.base_seed
    montecarlo_globals.ENABLE_RPACKET_TRACKING = (
        transport.enable_rpacket_tracking
    )
