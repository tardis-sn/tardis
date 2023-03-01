from tardis.io.atom_data.util import download_atom_data
from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation
from tardis.io.parsers.csvy import load_csvy
from astropy import units
from skopt.space import Real
from skopt.utils import use_named_args
from skopt import gp_minimize
import numpy as np
from functools import partial


CFG_FILE = "uniform_ejecta_model.yml"
TARGET_W = 0.5
NUM_ITERATIONS = 20


def get_simulation_error(cfg_file, v_inner):
    print(
        f"||||||||||||||||||||||||||||||||||||||||{v_inner}|||||||||||||||||||||||||||||||||||||||||||||||||||||"
    )
    cfg = Configuration.from_yaml(cfg_file)
    cfg.model.structure.velocity.start = (
        v_inner[0] * cfg.model.structure.velocity.start.unit
    )
    sim = Simulation.from_config(
        cfg, show_progress_bars=True, log_level="CRITICAL"
    )
    sim.run()
    return (TARGET_W - sim.plasma.w[0]) ** 2


if __name__ == "__main__":
    # Download the atomic data
    download_atom_data("kurucz_cd23_chianti_H_He")
    cfg = Configuration.from_yaml(CFG_FILE)

    START_VELOCITY = cfg.model.structure.velocity.start.value
    END_VELOCITY = cfg.model.structure.velocity.start.value + (
        (
            cfg.model.structure.velocity.stop.value
            - cfg.model.structure.velocity.start.value
        )
        * 0.9
    )

    # Search space
    SEARCH_SPACE = [
        Real(START_VELOCITY, END_VELOCITY, "uniform", name="v_inner")
    ]

    result = gp_minimize(
        partial(get_simulation_error, CFG_FILE),
        SEARCH_SPACE,
        n_calls=NUM_ITERATIONS,
    )

    print("######## FINAL RESULTS ##############")
    print(f"Minimum squared error: {result.fun}")
    # print(f"Best W: {BEST_W}")
    print(f"Best Inner Velocity: {result.x}")
