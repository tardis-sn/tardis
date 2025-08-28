"""
Progress bar utilities for Monte Carlo transport calculations.
"""

from __future__ import annotations

from tqdm.auto import tqdm

from tardis.util.environment import Environment

# Global progress bar instances - initialized at import time for thread safety
# Use tqdm.auto which automatically detects environment
iterations_pbar = tqdm(
    desc="Iterations:",
    bar_format="{desc:<}{bar}{n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]",
)
packet_pbar = tqdm(
    desc="Packets:   ",
    postfix="0",
    bar_format="{desc:<}{bar}{n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]",
)


def update_packet_pbar(
    i: int, current_iteration: int, no_of_packets: int, total_iterations: int
) -> None:
    """
    Update progress bars as each packet is propagated.

    Parameters
    ----------
    i : int
        Amount by which the progress bar needs to be updated.
    current_iteration : int
        Current iteration number.
    no_of_packets : int
        Total number of packets in one iteration.
    total_iterations : int
        Total number of iterations.
    """
    if packet_pbar.postfix == "":
        packet_pbar.postfix = "0"
    bar_iteration = int(packet_pbar.postfix) - 1

    # fix bar layout when run_tardis is called for the first time
    if iterations_pbar.total == None:
        fix_bar_layout(iterations_pbar, total_iterations=total_iterations)
    if packet_pbar.total == None:
        fix_bar_layout(packet_pbar, no_of_packets=no_of_packets)

    # display and reset progress bar when run_tardis is called again
    if iterations_pbar.n == total_iterations:
        if type(iterations_pbar).__name__ == "tqdm_notebook":
            iterations_pbar.container.close()
        fix_bar_layout(iterations_pbar, total_iterations=total_iterations)

    if bar_iteration > current_iteration:
        packet_pbar.postfix = current_iteration
        if type(packet_pbar).__name__ == "tqdm_notebook":
            # stop displaying last container
            packet_pbar.container.close()
        fix_bar_layout(packet_pbar, no_of_packets=no_of_packets)

    # reset progress bar with each iteration
    if bar_iteration < current_iteration:
        packet_pbar.reset(total=no_of_packets)
        packet_pbar.postfix = str(current_iteration + 1)

    packet_pbar.update(i)


def update_packets_pbar(i: int, no_of_packets: int) -> None:
    """
    Update packet progress bar with packet count updates only.

    This function only handles packet progress without iteration logic,
    suitable for use within the montecarlo main loop.

    Parameters
    ----------
    i : int
        Amount by which the progress bar needs to be updated.
    no_of_packets : int
        Total number of packets in one iteration.
    """
    # Initialize packet progress bar if needed
    if packet_pbar.total is None:
        fix_bar_layout(packet_pbar, no_of_packets=no_of_packets)

    packet_pbar.update(i)


def reset_packet_pbar(no_of_packets: int) -> None:
    """
    Reset the packet progress bar for a new iteration.

    Parameters
    ----------
    no_of_packets : int
        Total number of packets in the iteration.
    """
    packet_pbar.reset(total=no_of_packets)


def refresh_packet_pbar() -> None:
    """
    Refresh packet progress bar after each iteration.

    This function refreshes the visual display of the packet progress bar
    to ensure accurate rendering after each iteration completes.
    """
    packet_pbar.refresh()


def update_iterations_pbar(i: int) -> None:
    """
    Update progress bar for each iteration.

    Parameters
    ----------
    i : int
        Amount by which the progress bar needs to be updated.
    """
    iterations_pbar.update(i)


def initialize_iterations_pbar(total_iterations: int) -> None:
    """
    Initialize the iterations progress bar with the total number of iterations.

    Parameters
    ----------
    total_iterations : int
        Total number of iterations.
    """
    if iterations_pbar.total is None:
        fix_bar_layout(iterations_pbar, total_iterations=total_iterations)


def fix_bar_layout(
    bar, no_of_packets: int | None = None, total_iterations: int | None = None
) -> None:
    """
    Fix the layout of progress bars.

    This function handles the visual layout and configuration of progress bars,
    with proper environment detection for notebook vs terminal environments.

    Parameters
    ----------
    bar : tqdm.tqdm
        Progress bar instance to change the layout of.
    no_of_packets : int, optional
        Number of packets to be propagated, by default None.
    total_iterations : int, optional
        Total number of iterations, by default None.
    """
    # Check if we're in an environment that allows widget display and has notebook-style progress bar
    if (
        Environment.allows_widget_display()
        and type(bar).__name__ == "tqdm_notebook"
    ):
        bar.container = bar.status_printer(
            bar.fp,
            bar.total,
            bar.desc,
            bar.ncols,
        )
        if no_of_packets is not None:
            bar.reset(total=no_of_packets)
        if total_iterations is not None:
            bar.reset(total=total_iterations)

        # change the amount of space the prefix string of the bar takes
        # here, either packets or iterations
        bar.container.children[0].layout.width = "6%"

        # change the length of the bar
        bar.container.children[1].layout.width = "60%"

        # display the progress bar
        from IPython import display

        display.display(bar.container)
    else:
        if no_of_packets is not None:
            bar.reset(total=no_of_packets)
        if total_iterations is not None:
            bar.reset(total=total_iterations)
