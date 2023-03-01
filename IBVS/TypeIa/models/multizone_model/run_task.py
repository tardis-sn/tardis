from tardis.io.atom_data.util import download_atom_data
from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation
from tardis.io.parsers.csvy import load_csvy
from astropy import units
import numpy as np


def run_search(low_v, high_v, step_size, cfg_file):
    inner_velocities = []
    w_inner = []
    for v_inner in np.arange(low_v, high_v, step_size):
        inner_velocities.append(v_inner)

        cfg = Configuration.from_yaml(cfg_file)
        cfg.model.v_inner_boundary = v_inner * cfg.model.v_inner_boundary.unit

        sim = Simulation.from_config(
            cfg, show_progress_bars=True, log_level="CRITICAL"
        )
        sim.run()

        w_inner.append(sim.plasma.w[0])

        np.save("inner_velocities", np.array(inner_velocities))
        np.save("w_inner", np.array(w_inner))


if __name__ == "__main__":
    total_iterations = 0
    # Download the atomic data
    download_atom_data("kurucz_cd23_chianti_H_He")

    CFG_FILE = "multizone_model.yml"
    MODEL_FILE = "multizone_model.csvy"
    cfg = Configuration.from_yaml(CFG_FILE)
    csvy_model_config, csvy_model_data = load_csvy(MODEL_FILE)

    TARGET_W = 0.5

    START_VELOCITY = cfg.model.v_inner_boundary.value
    END_VELOCITY = (
        csvy_model_data.velocity.iloc[-1] * 1e5
    )  # Convert from km/s to cm/s

    # Global best results
    BEST_ERROR = np.inf
    BEST_W = -1
    BEST_VELOCITY = START_VELOCITY
    # Search space
    low_v = max(START_VELOCITY // 2, csvy_model_data.velocity.iloc[0] * 1e5)
    high_v = END_VELOCITY

    print("################ START OF SEARCH ##############################")

    while True:
        step_size = (high_v - low_v) // 5
        print(
            f"%%%%%%%%%%%%% Running search from {low_v} to {high_v} with steps {step_size}"
        )
        if step_size == 0:
            break

        run_search(low_v, high_v, step_size, CFG_FILE)

        inner_velocities = np.load("inner_velocities.npy")
        w_inner = np.load("w_inner.npy")

        # Set start and end for next level (Look around local best results)
        best_e = np.inf
        best_idx = -1
        for i, w in enumerate(w_inner):
            total_iterations += 1  # Increment iterations
            error = (w - TARGET_W) ** 2
            if error < BEST_ERROR:  # Update global solution
                BEST_ERROR = error
                BEST_W = w
                BEST_VELOCITY = inner_velocities[i]

            # Terminate early if error within range
            if BEST_ERROR < 1e-4:
                print("######## FINAL RESULTS ##############")
                print(f"Total Iterations: {total_iterations}")
                print(f"Minimum squared error: {BEST_ERROR}")
                print(f"Best W: {BEST_W}")
                print(f"Best Inner Velocity: {BEST_VELOCITY}")
                exit()

            if error < best_e:
                best_e = error
                best_idx = i

        if best_idx > 0:
            low_v = inner_velocities[best_idx - 1]

        if best_idx < len(inner_velocities) - 1:
            high_v = inner_velocities[best_idx + 1]

    print("######## FINAL RESULTS ##############")
    print(f"Total Iterations: {total_iterations}")
    print(f"Minimum squared error: {BEST_ERROR}")
    print(f"Best W: {BEST_W}")
    print(f"Best Inner Velocity: {BEST_VELOCITY}")
