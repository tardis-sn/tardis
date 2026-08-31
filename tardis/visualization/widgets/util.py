"""Utility classes and functions for widgets."""

import asyncio
import logging
import re
import panel as pn

logger = logging.getLogger(__name__)


class TableSummaryLabel:
    """
    Label like widget to show summary of a Panel table widget.

    Also handles aligning the label with the table columns exactly like a
    summary row.
    """

    def __init__(self, target_table, table_col_widths, label_key, label_value):
        """
        Initialize a TableSummaryLabel for a table widget.

        Parameters
        ----------
        target_table : PanelTableWidget
            Table widget whose summary label it is
        table_col_widths : list
            A list containing width of each column of table in order (including
            the index as 1st column). The width values must be proportions of
            100 i.e. they must sum to 100
        label_key : str
            Brief description of what label summarizes about table
        label_value : int
            Initial summary value of the table to be shown in label

        Notes
        -----
        TableSummaryLabel can only be created for a table with two
        columns (including index as 1st column) as of now
        """
        if len(table_col_widths) != 2:
            raise NotImplementedError(
                "Currently TableSummaryLabel can only be created for a table "
                "widget with exactly two columns (including index column)"
            )

        self.target_table = target_table
        self.table_col_widths = table_col_widths
        self.widget = self._create(label_key, label_value)

    def update_and_resize(self, value):
        """
        Update the label value and resize it as per the size of target table.

        Resizing is done in such a way so as to match the width of components
        of label with columns of target table, making it look like another row.
        This method should be called whenever there is any update in data or
        layout of target table.

        Parameters
        ----------
        value : int
            Value to be shown in label
        """
        # Update the value component with new value
        self.widget[1].object = f"<div style='padding:0px 2px; margin:0px;'>{str(value)}</div>"

        # For Panel components, sizing is handled automatically with sizing_mode='stretch_width'
        # No manual width calculation needed as Panel handles responsive sizing

    def get_value(self):
        """
        Get the current value displayed in the label.

        Returns
        -------
        str
            The current value as a string
        """
        # Extract value from HTML content
        html_content = self.widget[1].object
        match = re.search(r'<div[^>]*>([^<]*)</div>', html_content)
        if match:
            return match.group(1)
        return "0"

    def _create(self, key, value):
        """
        Create widget for TableSummaryLabel.

        Parameters
        ----------
        key : str
            Brief description of what label summarizes about target table
        value : int
            Initial summary value of the target table to be shown in label

        Returns
        -------
        panel.Row
            Widget containing all components of label
        """
        # Create Panel HTML components with styling
        key_component = pn.pane.HTML(
            f"<div style='text-align:right; padding:0px 2px; margin:0px;'> <b>{key}:</b> </div>",
            sizing_mode='stretch_width'
        )

        value_component = pn.pane.HTML(
            f"<div style='padding:0px 2px; margin:0px;'>{str(value)}</div>",
            sizing_mode='stretch_width'
        )

        return pn.Row(
            key_component,
            value_component,
            sizing_mode='stretch_width'
        )


class Timer:
    """Timer to implement debouncing using an asynchronous loop.

    Notes
    -----
    This class is reproduced from ipywidgets documentation, for more information
    please see https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20Events.html#debounce
    """

    def __init__(self, timeout, callback):
        """Initialize the Timer with delay time and delayed function.

        Parameters
        ----------
            timeout : float
            callback : function
        """
        self._timeout = timeout
        self._callback = callback

    async def _job(self):
        await asyncio.sleep(self._timeout)
        self._callback()

    def start(self):
        self._task = asyncio.ensure_future(self._job())

    def cancel(self):
        self._task.cancel()


def debounce(wait):
    """Decorator that will postpone a function's execution until after
     `wait` seconds have elapsed since the last time it was invoked.

    Parameters
    ----------
        wait : float

    Returns
    -------
        function

    Notes
    -----
    This decorator is reproduced from ipywidgets documentation, for more information
    please see https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20Events.html#debounce
    """

    def decorator(fn):
        timer = None

        def debounced(*args, **kwargs):
            nonlocal timer

            def call_it():
                fn(*args, **kwargs)

            if timer is not None:
                timer.cancel()
            timer = Timer(wait, call_it)
            timer.start()

        return debounced

    return decorator
