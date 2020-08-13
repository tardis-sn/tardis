import qgrid


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
