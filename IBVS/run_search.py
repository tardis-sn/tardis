from tardis.io.atom_data.util import download_atom_data
from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation
from tardis.io.atom_data.atom_web_download import get_atomic_repo_config
from tardis.io.config_internal import get_data_dir
from tardis import run_tardis
from astropy import units
import numpy as np
import argparse
import os


def run_simulation(cfg_file, v_inner, **sim_kwargs):
    """
    Runs a simulation with the required inner velocity

    Parameters
        ----------
        cfg_file : str Filepath of config file
        v_inner : float Desired inner velocity

        **sim_kwargs
            kwargs to pass to Simulation.from_config

        Returns
        -------
        tardis.montecarlo.Simulation
    """
    # Build config object
    cfg = Configuration.from_yaml(cfg_file)

    # Modify inner velocity
    if hasattr(cfg, "csvy_model"):
        cfg.model.v_inner_boundary = v_inner * cfg.model.v_inner_boundary.unit
    else:
        cfg.model.structure.velocity.start = (
            v_inner * cfg.model.structure.velocity.start.unit
        )

    # Run simulation
    atomic_data_name = get_atomic_repo_config()["default"]
    atom_data_path = os.path.join(get_data_dir(), f"{atomic_data_name}.h5")
    sim = run_tardis(cfg, atom_data=atom_data_path)
    sim = Simulation.from_config(cfg, **sim_kwargs)
    sim.run()

    return sim


def bayesian_optimization(
    cfg_file,
    v_inner_low,
    v_inner_high,
    target_w,
    error_threshold=None,
    max_iterations=None,
    prior="uniform",
    **sim_kwargs,
):
    """
    Runs Bayesian Optimization to search for the best v_inner

    Parameters
        ----------
        v_inner_low : float Lower bound of the search space
        v_inner_high : float Higher bound of the search space

        target_w : float Target value of innermost optical depth

        error_threshold : float Maximum allowed deviation of computed optical depth from target optical depth
        max_iterations : int Max number of iterations Bayesian Opt. should run for.
            If provided, Bayesian optimization stops if either this number is reached or error_threshold is obeyed,
                whichever is earlier
            If not provided, Bayesian optimization keeps on running till error_threshold is obeyed
        prior : "uniform" or "log-uniform", default="uniform" , Prior for sampling from search space

        **sim_kwargs
            kwargs to pass to Simulation.from_config

        Returns
        -------
        (float, tardis.montecarlo.Simulation) : (Best inner velocity, simulation run on this best velocity)
    """
    if error_threshold is None and max_iterations is None:
        raise ValueError(
            "Error threshold and max iterations both cannot be None"
        )

    from skopt import gp_minimize
    from skopt.callbacks import DeltaYStopper
    from skopt.space import Real

    # Declare search space
    search_space = [
        Real(v_inner_low, v_inner_high, prior=prior, name="v_inner")
    ]

    # Early stopping callback
    callback = (
        DeltaYStopper(error_threshold, n_best=1)
        if error_threshold is not None
        else error_threshold
    )

    def compute_error(v_inner):
        """Computes objective for Bayesian optimization to minimize"""
        sim = run_simulation(cfg_file, v_inner[0], **sim_kwargs)
        return abs(target_w - sim.plasma.w[0])

    result = gp_minimize(
        compute_error,
        search_space,
        callback=callback,
        n_calls=max_iterations,
        verbose=True,
        n_jobs=-1,
    )

    # Run simulation with the best found v_inner
    best_v_inner = result.x
    best_sim = run_simulation(cfg_file, best_v_inner, **sim_kwargs)

    return best_v_inner, best_sim


def zoom_in_search(
    cfg_file,
    v_inner_low,
    v_inner_high,
    target_w,
    num_dividers=5,
    error_threshold=None,
    max_levels=None,
    **sim_kwargs,
):
    """
    Zoom-in Search to search for the best v_inner
    The function runs in a recursive fashion

    Parameters
        ----------
        v_inner_low : float Lower bound of the search space
        v_inner_high : float Higher bound of the search space

        target_w : float Target value of innermost optical depth

        num_dividers : int Number of equally-spaced point to search at each level
        error_threshold : float Maximum allowed deviation of computed optical depth from target optical depth
        max_levels : int Max number of levels to keep zooming-in
            If provided, the algorithm stops if either this level is reached or error_threshold is obeyed,
                whichever is earlier
            If not provided, the algorithm keeps on running till error_threshold is obeyed

        **sim_kwargs
            kwargs to pass to Simulation.from_config

        Returns
        -------
        (float, tardis.montecarlo.Simulation) : (Best inner velocity, simulation run on this best velocity)
    """
    # Global best estimates
    best_error = np.inf
    best_v_inner = -1
    best_sim = None

    # Runs simulations for all dividers in this level
    v_inner_candidates = np.linspace(v_inner_low, v_inner_high, num_dividers)
    best_idx = -1
    for idx, v_inner in enumerate(v_inner_candidates):
        sim = run_simulation(cfg_file, v_inner, **sim_kwargs)
        w = sim.plasma.w[0]

        # Find index with minimum error
        if abs(w - target_w) < best_error:
            best_error = abs(w - target_w)
            best_idx = idx
            best_sim = sim

    best_v_inner = v_inner_candidates[best_idx]

    # Base case: Return best result from this level is base case criteria matches
    if (error_threshold is not None and best_error <= error_threshold) or (
        max_levels is not None and max_levels == 0
    ):
        return best_v_inner, best_sim

    # Recurse to lower levels
    v_inner_low_next = (
        v_inner_candidates[best_idx - 1]
        if best_idx > 0
        else v_inner_candidates[0]
    )
    v_inner_high_next = (
        v_inner_candidates[best_idx + 1]
        if best_idx < len(v_inner_candidates) - 1
        else v_inner_candidates[-1]
    )
    best_v_inner_next, best_sim_next = zoom_in_search(
        cfg_file,
        v_inner_low_next,
        v_inner_high_next,
        target_w,
        num_dividers=num_dividers,
        error_threshold=error_threshold,
        max_levels=max_levels - 1
        if max_levels is not None
        else None,  # Decrease max levels by 1
        **sim_kwargs,
    )
    best_error_next = abs(best_sim_next.plasma.w[0] - target_w)

    # Update best numbers if lower levels had them better
    if best_error_next < best_error:
        best_error = best_error_next
        best_sim = best_sim_next
        best_v_inner = best_v_inner_next

    return best_v_inner, best_sim


def run_search(low_v, high_v, step_size, cfg_file):
    inner_velocities = []
    w_inner = []
    for v_inner in np.arange(low_v, high_v, step_size):
        inner_velocities.append(v_inner)

        cfg = Configuration.from_yaml(cfg_file)
        cfg.model.structure.velocity.start = (
            v_inner * cfg.model.structure.velocity.start.unit
        )

        sim = Simulation.from_config(
            cfg, show_progress_bars=True, log_level="CRITICAL"
        )
        sim.run()

        w_inner.append(sim.plasma.w[0])

        np.save("inner_velocities", np.array(inner_velocities))
        np.save("w_inner", np.array(w_inner))


if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(
        prog="run_v_inner_search",
        description="Searches v_inner which gives a target optical depth of innermost shell",
    )

    parser.add_argument(
        "--atom-data", type=str, default="kurucz_cd23_chianti_H_He"
    )
    parser.add_argument(
        "--config",
        "--cfg",
        type=str,
        default="/Users/ansh/Documents/tardis/IBVS/TypeIa/models/multizone_model/multizone_model.yml",
    )
    parser.add_argument(
        "--target_w", type=float, default=0.5, help="Target optical depth"
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.01,
        help="Maximum error of target_w allowed",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=40,
        help="Max number of iterations for Bayesian Opt",
    )
    parser.add_argument(
        "--levels", type=int, help="Max depth of zoom-in search"
    )
    parser.add_argument(
        "--algorithm",
        "--algo",
        type=str,
        choices=["zoom", "bayesian"],
        default="bayesian",
    )
    parser.add_argument(
        "--dividers",
        type=int,
        default=5,
        help="Number of dividers for zoom-in search",
    )
    parser.add_argument(
        "--prior",
        type=str,
        choices=["uniform", "log-uniform"],
        default="uniform",
        help="Prior for sampling from search space for Bayesian Opt",
    )

    args = parser.parse_args()

    # Download the atomic data
    download_atom_data(args.atom_data)

    cfg = Configuration.from_yaml(args.config)

    if hasattr(cfg, "csvy_model"):
        from tardis.io.parsers.csvy import load_csvy

        # Get original inner velocity
        start_velocity = cfg.model.v_inner_boundary
        # Get outer velocity
        csvy_model_config, csvy_model_data = load_csvy(cfg.csvy_model)
        end_velocity_value = csvy_model_data.velocity.iloc[-1]
        end_velocity_unit = units.Unit(
            list(
                filter(
                    lambda row: row.get("name") == "velocity",
                    csvy_model_config.get("datatype").get("fields"),
                )
            )[0].get("unit")
        )
        # Convert end velocity to same unit as start velocity
        end_velocity = (end_velocity_value * end_velocity_unit).to(
            start_velocity.unit
        )
    else:
        start_velocity = cfg.model.structure.velocity.start
        end_velocity = cfg.model.structure.velocity.stop

    # Search space
    low_v = start_velocity.value  # TODO: Set correct lower bound
    high_v = end_velocity.value

    # Run search
    if args.algorithm == "zoom":
        best_v_inner, best_sim = zoom_in_search(
            args.config,
            low_v,
            high_v,
            args.target_w,
            num_dividers=args.dividers,
            error_threshold=args.threshold,
            max_levels=args.levels,
            show_progress_bars=True,
        )
    elif args.algorithm == "bayesian":
        best_v_inner, best_sim = bayesian_optimization(
            args.config,
            low_v,
            high_v,
            args.target_w,
            error_threshold=args.threshold,
            max_iterations=args.iterations,
            prior=args.prior,
            show_progress_bars=True,
        )
    else:
        raise ValueError("Wrong Algorithm")

    print("######## FINAL RESULTS ##############")
    print(f"Best W: {best_sim.plasma.w[0]}")
    print(f"Best Inner Velocity: {best_v_inner}")
