"""Utility classes and functions for widgets."""

import asyncio
import logging

import ipywidgets as ipw
import pandas as pd
import panel as pn

logger = logging.getLogger(__name__)

pn.extension('tabulator')

def create_table_widget(
    data, col_widths, table_options=None, changeable_col=None
):
    """
    Create table widget object which supports interaction and updating the data.

    Parameters
    ----------
    data : pandas.DataFrame
        Data you want to display in table widget
    col_widths : list
        A list containing width of each column of data in order (including
        the index as 1st column). The width values must be proportions of
        100 i.e. they must sum to 100.
    table_options : dict, optional
        A dictionary to specify options to use when creating interactive table
        widget (same as :grid_options: of qgrid, specified in Notes section of
        their `API documentation <https://qgrid.readthedocs.io/en/latest/#qgrid.show_grid>_`).
        Invalid keys will have no effect and valid keys will override or add to
        the options used by this function.
    changeable_col : dict, optional
        A dictionary to specify the information about column which will
        change its name when data in generated table widget updates. It
        must have two keys - :code:`index` to specify index of changeable
        column in dataframe :code:`data` as an integer, and :code:`other_names`
        to specify all possible names changeable column will get as a list
        of strings. Default value :code:`None` indicates that there is no
        changable column.

    Returns
    -------
    panel.widgets.Tabulator
        Table widget object
    """

    # Check whether passed col_widths list is correct or not
    if len(col_widths) != data.shape[1] + 1:
        raise ValueError(
            "Size of column widths list do not match with "
            "number of columns + 1 (index) in dataframe"
        )

    # Note: Since forceFitColumns is enabled by default in grid_options,
    # the column widths (when all specified) get applied in proportions,
    # despite their original unit is px thus it's better they sum to 100
    if sum(col_widths) != 100:
        raise ValueError(
            "Column widths are not proportions of 100 (i.e. "
            "they do not sum to 100)"
        )
    
        # Convert qgrid-style options to Tabulator options
    tabulator_options = {
        "sortable": False,
        "selectable": False,
        "show_index": True,
    }
    
    # Handle pagination options
    if table_options:
        if "maxVisibleRows" in table_options:
            tabulator_options["pagination"] = "local"
            tabulator_options["page_size"] = table_options["maxVisibleRows"]
        if "editable" in table_options:
            tabulator_options["selectable"] = table_options["editable"]
        if "sortable" in table_options:
            tabulator_options["sortable"] = table_options["sortable"]
        if "filterable" in table_options:
            # Panel's Tabulator doesn't have a direct filterable option
            # But we can add filters to each column
            pass

    # Use default width for better responsiveness (units in pixels)
    default_width = 800

    # Preparing dictionary that defines column widths
    column_widths = {}
    cols_with_index = [data.index.name] + data.columns.to_list()
    for col_name, width_prop in zip(cols_with_index, col_widths):
        column_widths[col_name] = int(width_prop / 100 * default_width)

    # We also need to define widths for different names of changeable column
    # Handle changeable column if specified
    if changeable_col and {"index", "other_names"}.issubset(set(changeable_col.keys())):
        changeable_width = col_widths[changeable_col["index"]]
        for col_name in changeable_col["other_names"]:
            column_widths[col_name] = int(changeable_width / 100 * default_width)

    # Create the table widget using Panel's Tabulator
    tabulator = pn.widgets.Tabulator(
        data,
        widths=column_widths,
        width=default_width,
        **tabulator_options
    )
    
    return tabulator



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
        target_table : qgrid.QgridWidget
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


        Notes
        -----
        The width resizing operation is adapted for Panel Tabulator tables.
        """
        self.widget.children[1].value = str(value)

        try:
            # Get the width of the Panel Tabulator widget
            table_width = self.target_table.width
            
            # Convert to integer if it's a string with 'px'
            if isinstance(table_width, str) and table_width.endswith('px'):
                table_width = int(table_width.rstrip('px'))
            elif not isinstance(table_width, (int, float)):
                # Use a default width if not defined
                table_width = 800
                
        except AttributeError:
            logger.warning(
                "target_table doesn't have a width attribute, label "
                "cannot be resized!",
                exc_info=1,
            )
            return

        # Check if pagination is enabled
        has_pagination = getattr(self.target_table, 'pagination', None) == 'local'
        pagination_height = 30 if has_pagination else 0
        
        # Panel's Tabulator handles scrollbars differently than qgrid
        # No need to subtract scrollbar width explicitly
        
        # Distribute the table width in proportions of column width to label components
        self.widget.children[0].layout.width = f"{(table_width) * self.table_col_widths[0]/100}px"
        self.widget.children[1].layout.width = f"{(table_width) * self.table_col_widths[1]/100}px"

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
        ipywidgets.Box
            Widget containing all componets of label
        """
        # WARNING: Use dictionary instead of ipw.Layout for specifying layout
        # of ipywidgets, otherwise there will be unintended behavior
        component_layout_options = dict(
            flex="0 0 auto",  # to prevent shrinking of flex-items
            padding="0px 2px",  # match with header
            margin="0px",  # remove default 1px margin
        )

        return ipw.Box(
            [
                ipw.HTML(
                    f"<div style='text-align:right;'> <b>{key}:<b> </div>",
                    layout=component_layout_options,
                ),
                ipw.HTML(
                    str(value),
                    layout=component_layout_options,
                ),
            ],
            layout=dict(
                display="flex",
                justify_content="flex-start",
            ),
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
