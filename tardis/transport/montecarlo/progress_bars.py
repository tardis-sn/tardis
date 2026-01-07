"""
Progress bar utilities for Monte Carlo transport calculations.
"""

from tqdm.auto import tqdm

from tardis.util.environment import Environment

# Global progress bar instances - lazily initialized to avoid output at import time
_iterations_pbar = None
_packet_pbar = None


def _get_iterations_pbar():
    """
    Lazily initialize and return the iterations progress bar.

    Returns
    -------
    tqdm.tqdm
        The iterations progress bar instance.
    """
    global _iterations_pbar
    if _iterations_pbar is None:
        _iterations_pbar = tqdm(
            desc="Iterations:",
            bar_format="{desc:<}{bar}{n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]",
        )
    return _iterations_pbar


def _get_packet_pbar():
    """
    Lazily initialize and return the packet progress bar.

    Returns
    -------
    tqdm.tqdm
        The packet progress bar instance.
    """
    global _packet_pbar
    if _packet_pbar is None:
        _packet_pbar = tqdm(
            desc="Packets:   ",
            postfix="0",
            bar_format="{desc:<}{bar}{n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]",
        )
    return _packet_pbar


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
    packet_pbar = _get_packet_pbar()
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
    _get_packet_pbar().reset(total=no_of_packets)


def refresh_packet_pbar() -> None:
    """
    Refresh packet progress bar after each iteration.

    This function refreshes the visual display of the packet progress bar
    to ensure accurate rendering after each iteration completes.
    """
    _get_packet_pbar().refresh()


def update_iterations_pbar(i: int) -> None:
    """
    Update progress bar for each iteration.

    Parameters
    ----------
    i : int
        Amount by which the progress bar needs to be updated.
    """
    _get_iterations_pbar().update(i)


def initialize_iterations_pbar(total_iterations: int) -> None:
    """
    Initialize the iterations progress bar with the total number of iterations.

    Parameters
    ----------
    total_iterations : int
        Total number of iterations.
    """
    iterations_pbar = _get_iterations_pbar()
    if iterations_pbar.total is None:
        fix_bar_layout(iterations_pbar, total_iterations=total_iterations)


def close_progress_bars():
    """
    Close all progress bars if they have been initialized.

    This function is used for cleanup, e.g., at the end of test sessions.
    """
    global _iterations_pbar, _packet_pbar
    if _packet_pbar is not None:
        _packet_pbar.close()
    if _iterations_pbar is not None:
        _iterations_pbar.close()


def fix_bar_layout(bar, no_of_packets=None, total_iterations=None):
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
