import logging
import qgrid
import ipywidgets as ipw

logger = logging.getLogger(__name__)


def create_table_widget(
    data, col_widths, table_options=None, changeable_col=None
):
    """Create table widget object which supports interaction and updating
    the data

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
    qgrid.QgridWidget
        Table widget object
    """
    # Setting the options to be used for creating table widgets
    grid_options = {
        "sortable": False,
        "filterable": False,
        "editable": False,
        "minVisibleRows": 2,
    }
    if table_options:
        grid_options.update(table_options)

    column_options = {
        "minWidth": None,
    }

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

    # Preparing dictionary that defines column widths
    cols_with_index = [data.index.name] + data.columns.to_list()
    column_widths_definitions = {
        col_name: {"width": col_width}
        for col_name, col_width in zip(cols_with_index, col_widths)
    }

    # We also need to define widths for different names of changeable column
    if changeable_col:
        if {"index", "other_names"}.issubset(set(changeable_col.keys())):
            column_widths_definitions.update(
                {
                    col_name: {"width": col_widths[changeable_col["index"]]}
                    for col_name in changeable_col["other_names"]
                }
            )
        else:
            raise ValueError(
                "Changeable column dictionary does not contain "
                "'index' or 'other_names' key"
            )

    # Create the table widget using qgrid
    return qgrid.show_grid(
        data,
        grid_options=grid_options,
        column_options=column_options,
        column_definitions=column_widths_definitions,
    )


class TableSummaryLabel:
    """Label like widget to show summary of a qgrid table widget. Also handles
    aligning the label with the table columns exactly like a summary row.
    """

    def __init__(self, target_table, table_col_widths, label_key, label_value):
        if len(table_col_widths) != 2:
            raise NotImplementedError(
                "Currently TableSummaryLabel can only be created for a table "
                "widget with exactly two columns (including index column)"
            )

        self.target_table = target_table
        self.table_col_widths = table_col_widths
        self.widget = self._create(label_key, label_value)

    def update_and_resize(self, value):
        """Update the label value and resize the width of label components to
        align them with the columns of target table. This should be called
        whenever there is any update in target table.

        Parameters
        ----------
        value : int
            Value to be shown in label

        Note
        ----
        The width resizing operation is highly dependent on qgrid tables' 
        layout. So it may not remain accurate if there happens any CSS change 
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

        if (
            len(self.target_table.df)
            > self.target_table.grid_options["maxVisibleRows"]
        ):
            table_width -= 12  # 12px is space consumed by scroll bar of qgrid

        self.widget.children[
            0
        ].layout.width = f"{(table_width) * self.table_col_widths[0]/100}px"
        self.widget.children[
            1
        ].layout.width = f"{(table_width) * self.table_col_widths[1]/100}px"

    def _create(self, key, value):
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
                ipw.HTML(str(value), layout=component_layout_options,),
            ],
            layout=dict(display="flex", justify_content="flex-start",),
        )
