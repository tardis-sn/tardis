try:
    import lineid_plot
    from lineid_plot import get_line_flux

except ImportError:
    print("Please pip install lineid_plot to use the lineid_plotter function")
    raise


def lineid_plotter(
    ax,
    line_wavelengths,
    line_labels,
    spectrum_wavelengths,
    spectrum_data,
    style="top",
    plotter_kwargs=None,
    lineid_kwargs=None,
):
    """Use the lineid_plot package to plot user specified lines on a plotter's spectrum.

    Parameters
    ----------
    ax : matplotlib.axes._subplots.AxesSubplot
        Axes for which the spectral lines will be plotted.
    line_wavelengths : list of float
        Wavelength values of the lines to be plotted
    line_labels : list of str
        Labels to be used for the lines at the wavelength positions
    spectrum_wavelengths : astropy.units.Quantity
        Wavelengths of the spectrum where the lines will be plotted.
    spectrum_data : astropy.units.Quantity
        Fluxes or luminosity of the spectrum where the lines will be plotted.
    style : {'top', 'inside', 'along'}, optional
        Preset styles for the lines markers. (default: 'top')
        If 'top', the lines will be plotted above the top axes.
        If 'inside', the lines will be plotted inside the axes.
        If 'along', the lines will be plotted along the spectrum.
    lineid_kwargs : dict, optional
        kwargs passed to the lineid_plot.plot_line_ids method, by default None

    Returns
    -------
    fig, ax: matplotlib.figure.Figure, matplotlib.axes._subplots.AxesSubplot
        figure and axes with the spectral lines marked.
    """

    if plotter_kwargs is None:
        plotter_kwargs = {}

    if lineid_kwargs is None:
        lineid_kwargs = {}

    # if len(line_wavelengths) != len(line_labels):
    #     raise ValueError(
    #         "line_wavelengths and line_labels must have the same length"
    #     )

    if style == "top":
        lineid_plot.plot_line_ids(
            spectrum_wavelengths, spectrum_data, line_wavelengths, line_labels, ax=ax, **lineid_kwargs
        )

    elif style == "inside":
        box_loc = ax.transData.inverted().transform(
            ax.transAxes.transform((0, 0.9))
        )[1]
        arrow_loc = ax.transData.inverted().transform(
            ax.transAxes.transform((0, 0.8))
        )[1]

        lineid_plot.plot_line_ids(
            spectrum_wavelengths,
            spectrum_data,
            line_wavelengths,
            line_labels,
            ax=ax,
            arrow_tip=arrow_loc,
            box_loc=box_loc,
            **lineid_kwargs,
        )

    elif style == "along":

        # in data the center of boxs
        box_loc = ax.transData.inverted().transform(
            ax.transAxes.transform((0, 0.9))
        )[1]
        arrow_len = (
            box_loc
            - ax.transData.inverted().transform(
                ax.transAxes.transform((0, 0.8))
            )[1]
        )
        # 1/10 of the axes length, but measured up here because 0 to 0.1 in axes in the negative luminosity
        arrow_loc = box_loc - arrow_len  # in data, the start of the arrow
        marker_len = 2 * arrow_len
        fluxes = get_line_flux(
            line_wavelengths, spectrum_wavelengths.value[::-1], spectrum_data.value[::-1]
        )  # uses np.interp needs to monotonically inc

        # fix the max height to be that of the inside style
        arrow_tips = fluxes + marker_len
        arrow_tips[arrow_tips > 0.8] = 0.8
        box_locs = fluxes + marker_len + arrow_len
        box_locs[box_locs > 0.9] = 0.9

        lineid_plot.plot_line_ids(
            spectrum_wavelengths,
            spectrum_data,
            line_wavelengths,
            line_labels,
            ax=ax,
            arrow_tip=list(arrow_tips),
            box_loc=list(box_locs),
            **lineid_kwargs,
        )
    else:
        raise ValueError(
            "style must be one of 'top', 'inside', or 'along'"
        )

    return ax
