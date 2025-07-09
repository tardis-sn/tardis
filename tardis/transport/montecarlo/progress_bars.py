"""
Progress bar utilities for Monte Carlo transport calculations.
"""

from typing import Optional

from tqdm.auto import tqdm

# Global progress bar instances
iterations_pbar: Optional[tqdm] = None
packet_pbar: Optional[tqdm] = None


def _initialize_progress_bars() -> None:
    """
    Initialize progress bars if they haven't been created yet.

    This function lazily initializes the global progress bar instances,
    creating them only when first needed to avoid early tqdm initialization.
    """
    global iterations_pbar, packet_pbar

    if iterations_pbar is None or packet_pbar is None:
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
    # Initialize progress bars if needed
    _initialize_progress_bars()

    # These assertions help the type checker understand these are not None after initialization
    assert iterations_pbar is not None
    assert packet_pbar is not None

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


def refresh_packet_pbar() -> None:
    """
    Refresh packet progress bar after each iteration.

    This function refreshes the visual display of the packet progress bar
    to ensure accurate rendering after each iteration completes.
    """
    # Initialize progress bars if needed
    _initialize_progress_bars()
    assert packet_pbar is not None
    packet_pbar.refresh()


def update_iterations_pbar(i: int) -> None:
    """
    Update progress bar for each iteration.

    Parameters
    ----------
    i : int
        Amount by which the progress bar needs to be updated.
    """
    # Initialize progress bars if needed
    _initialize_progress_bars()
    assert iterations_pbar is not None
    iterations_pbar.update(i)


def fix_bar_layout(
    bar: tqdm,
    no_of_packets: Optional[int] = None,
    total_iterations: Optional[int] = None,
) -> None:
    """
    Fix the layout of progress bars.

    This function handles the visual layout and configuration of progress bars,
    with automatic detection of notebook vs terminal environments via tqdm.auto.

    Parameters
    ----------
    bar : tqdm.tqdm
        Progress bar instance to change the layout of.
    no_of_packets : int, optional
        Number of packets to be propagated, by default None.
    total_iterations : int, optional
        Total number of iterations, by default None.
    """
    if no_of_packets is not None:
        bar.reset(total=no_of_packets)
    if total_iterations is not None:
        bar.reset(total=total_iterations)

    # Handle notebook-specific layout if applicable
    if hasattr(bar, "container") and hasattr(bar.container, "children"):
        try:
            # change the amount of space the prefix string of the bar takes
            bar.container.children[0].layout.width = "6%"
            # change the length of the bar
            bar.container.children[1].layout.width = "60%"

            # Note: tqdm.auto automatically handles display, no need to manually display
        except (AttributeError, ImportError):
            # Not in a notebook environment or layout not available
            pass
