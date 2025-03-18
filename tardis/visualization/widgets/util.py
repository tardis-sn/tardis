"""Utility classes and functions for widgets."""

import asyncio
import logging
import panel as pn

logger = logging.getLogger(__name__)

pn.extension('tabulator')

def create_table_widget(
    data, col_widths, table_options=None, changeable_col=None
):
    """
    Create an interactive table widget using Panel's Tabulator, supporting
    data display, interaction, and optional column name changes.

    Parameters
    ----------
    data : pandas.DataFrame
        Data you want to display in table widget
    col_widths : list
        A list containing width of each column of data in order (including
        the index as 1st column). The width values must be proportions of
        100 i.e. they must sum to 100.
    table_options : dict, optional
        A dictionary specifying configuration options for the Tabulator widget.
        Overrides the default options where applicable.
    changeable_col : dict, optional
        A dictionary specifying the information about column that may change its name when
        the data updates. It must contain two keys:
        - :code:`index`: The index of the column in the DataFrame :code:`data`.
        - :code:`other_names`: A list of possible new names for the column.

    Returns
    -------
    panel.widgets.Tabulator
        An interactive Tabulator table widget displaying the DataFrame.

    Raises
    ------
    ValueError
        If the length of :code:`col_widths` does not match the number of
        columns + 1 (for the index), or if the column widths do not sum to 100.
    ValueError
        If :code:`changeable_col` does not contain both 'index' and 'other_names' keys.
    """
    if len(col_widths) != data.shape[1] + 1:
        raise ValueError(
            "Size of column widths list do not match with "
            "number of columns + 1 (index) in dataframe"
        )
    if sum(col_widths) != 100:
        raise ValueError(
            "Column widths are not proportions of 100 (i.e. "
            "they do not sum to 100)"
        )

    # Default Tabulator options
    tabulator_options = {
        'layout': 'fit_data_fill',
        'pagination': None,
        'selectable': 1,
    }
    num_rows = data.shape[0]
    if num_rows > 20:
        tabulator_options['height'] = 550
    else:
        tabulator_options['height'] = None
    if table_options:
        tabulator_options.update(table_options)

    # Define widths for columns (including index)
    widths = {data.index.name or 'index': f'{col_widths[0]}%'}  # Handle case where index.name is None
    widths.update({col: f'{col_widths[i+1]}%' for i, col in enumerate(data.columns)})
    custom_css = """
    .tabulator-header {
        height: auto !important;
        min-height: 40px;  /* Ensure enough height for headers */
        white-space: normal;  /* Allow wrapping if needed */
        overflow: visible !important;  /* Prevent clipping */
        text-overflow: clip;  /* Prevent truncation */
    }
    .tabulator-col-title {
        white-space: normal;  /* Allow wrapping */
        overflow: visible !important;
        text-overflow: clip;
        padding: 4px;  /* Add padding for readability */
    }
    .tabulator-tableholder {
        overflow-y: auto !important;  /* Ensure scrollbar works */
    }
    """

    if changeable_col:
        if not {"index", "other_names"}.issubset(set(changeable_col.keys())):
            raise ValueError(
                "Changeable column dictionary does not contain "
                "'index' or 'other_names' key"
            )

    return pn.widgets.Tabulator(
        data,
        **tabulator_options,
        widths=widths,
        stylesheets=[
            ":host {--mdc-ripple-color: transparent;}",
            custom_css
        ]
    )


class TableSummaryLabel:
    """
    Label like widget to show summary of a qgrid table widget.

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
        The width resizing operation is highly dependent on qgrid tables'
        layout. So it may not remain precise if there happens any CSS change
        in upcoming versions of qgrid.
        """
        self.widget.children[1].value = str(value)

        try:
            table_width = int(self.target_table.layout.width.rstrip("px"))
        except AttributeError:
            logger.warning(
                "target_table doesn't have any fixed width defined, label "
                "cannot be resized!",
                exc_info=1,
            )
            return

        max_rows_allowed = self.target_table.grid_options["maxVisibleRows"]
        if len(self.target_table.df) > max_rows_allowed:
            table_width -= 12  # 12px is space consumed by scroll bar of qgrid

        # Distribute the table width in proportions of column width to label components
        self.widget.children[
            0
        ].layout.width = f"{(table_width) * self.table_col_widths[0]/100}px"
        self.widget.children[
            1
        ].layout.width = f"{(table_width) * self.table_col_widths[1]/100}px"

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
