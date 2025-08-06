"""Utility classes and functions for widgets."""

import asyncio
import logging
import re
import panel as pn

logger = logging.getLogger(__name__)


class PanelTableWidget:
    """Panel-based table widget that mimics qgrid functionality."""

    def __init__(self, data, table_options=None):
        """Initialize the Panel table widget.

        Parameters
        ----------
        data : pandas.DataFrame
            Data to display in table widget
        table_options : dict, optional
            Table configuration options
        """
        self._df = data.copy()
        self.table_options = table_options or {}
        self._selected_rows = []
        self._selection_callbacks = []

        # Create the Panel table
        self._create_table()

    def _create_table(self):
        """Create the Panel table widget."""
        # Configure table parameters
        pagination = 'remote' if self.table_options.get('maxVisibleRows') else None
        page_size = self.table_options.get('maxVisibleRows', len(self._df))

        # Create the table
        self.table = pn.widgets.Tabulator(
            self._df,
            pagination=pagination,
            page_size=page_size,
            selectable=True,  # Single row selection (radio button style)
            show_index=True,
            sizing_mode='stretch_width',
            height=min(400, max(200, len(self._df) * 30 + 50)),
            disabled=True,  # Make cells non-editable
            configuration={
                'selectable': 'highlight',  # Make entire row selectable/highlightable
                'selectableRangeMode': 'click',  # Allow clicking anywhere on row to select
            }
        )

        # Set up selection callback
        self.table.param.watch(self._on_selection_change, 'selection')

    def _on_selection_change(self, event):
        """Handle selection changes in the table."""
        if event.new:
            # With single selection, Panel sends a list with one item or just the index
            if isinstance(event.new, list):
                selected_indices = [int(i) for i in event.new]
            else:
                selected_indices = [int(event.new)]
            self._selected_rows = selected_indices
        else:
            self._selected_rows = []
            selected_indices = []

        # Call registered callbacks
        for callback in self._selection_callbacks:
            # Create event dict similar to qgrid format
            event_dict = {
                'new': selected_indices,
                'old': [int(i) for i in getattr(event, 'old', [])] if hasattr(event, 'old') and event.old else [],
                'source': 'user'
            }
            callback(event_dict, self)

    @property
    def df(self):
        """Get the current dataframe."""
        return self._df

    @df.setter
    def df(self, value):
        """Set the dataframe and update the table."""
        self._df = value.copy()
        self.table.value = self._df

    def get_selected_rows(self):
        """Get the currently selected row indices."""
        return self._selected_rows

    def change_selection(self, index_values):
        """Programmatically change the selection by index values."""
        if not index_values:
            self.table.selection = []
            self._selected_rows = []
        else:
            # For single selection, only take the first index value
            idx_val = index_values[0] if isinstance(index_values, list) else index_values
            try:
                row_pos = self._df.index.get_loc(idx_val)
                row_position = int(row_pos)  # Convert to Python int
                self.table.selection = [row_position]  # Panel expects a list even for single selection
                self._selected_rows = [row_position]
            except (KeyError, TypeError) as e:
                print(f"Error selecting index {idx_val}: {e}")
                self.table.selection = []
                self._selected_rows = []

    def on(self, event_type, callback):
        """Register an event callback."""
        if event_type == 'selection_changed':
            self._selection_callbacks.append(callback)

    @property
    def layout(self):
        """Get the layout object for compatibility."""
        return self.table.param

    def _repr_mimebundle_(self, include=None, exclude=None):
        """Display the table widget."""
        return self.table._repr_mimebundle_(include, exclude)

    def show(self):
        """Show the table widget using Panel's show method."""
        return self.table.show()


def create_table_widget(data, table_options=None):
    """
    Create table widget object which supports interaction and updating the data.

    Parameters
    ----------
    data : pandas.DataFrame
        Data you want to display in table widget
    table_options : dict, optional
        A dictionary to specify options to use when creating interactive table
        widget. Supported options: maxVisibleRows.

    Returns
    -------
    PanelTableWidget
        Table widget object
    """
    return PanelTableWidget(data, table_options)


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
