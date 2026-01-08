"""
Tests for progress bar utilities.
"""

import tardis.transport.montecarlo.progress_bars as pb


def test_close_progress_bars_after_initialization():
    """Test that close_progress_bars properly closes initialized progress bars."""
    # Initialize the progress bars by calling the getter functions
    pb._get_iterations_pbar()
    pb._get_packet_pbar()

    # Verify they are initialized
    assert pb._iterations_pbar is not None
    assert pb._packet_pbar is not None

    # Close them
    pb.close_progress_bars()

    # Reset globals for other tests
    pb._iterations_pbar = None
    pb._packet_pbar = None
