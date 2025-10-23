import logging

from tardis.io.logger.logger import logging_state


def test_logging_state_without_debug_section():
    # Simulate a minimal configuration object without a debug section.
    # The function should not raise and should return a logger instance.
    config = {}
    cols, tardislogger = logging_state(None, config, specific_log_level=None, display_logging_widget=False)
    assert tardislogger is not None
    # If widgets unavailable, cols may be {} or None depending on environment
    assert (cols is None) or isinstance(cols, dict)
