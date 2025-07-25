import logging
import re

import pandas as pd


class PanelLoggingHandler(logging.Handler):
    def __init__(self, feeds):
        super().__init__()
        self.feeds = feeds
        self.feeds[0].append("Logger init")

    @staticmethod
    def _remove_ansi_escape_sequences(text):
        """Remove ANSI escape sequences from string.

        Parameters
        ----------
        text : str
            The text containing ANSI escape sequences.

        Returns
        -------
        str
            Cleaned text with ANSI escape sequences removed.
        """
        ansi_escape = re.compile(r"\x1B[@-_][0-?]*[ -/]*[@-~]")
        return ansi_escape.sub("", text)

    def emit(self, record):
        level = int(record.levelno / 10)
        if isinstance(record.msg, pd.DataFrame):
            record = record.msg
        else:
            record = self._remove_ansi_escape_sequences(self.format(record))
        self.feeds[level].append(record)
        self.feeds[0].append(record)
