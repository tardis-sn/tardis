try:
    import lineid_plot
    from lineid_plot import get_line_flux

except ImportError:
    print("Please pip install lineid_plot to use the lineid_plotter function")
    raise


def lineid_plotter(
    plotter,
    wavelengths,
    line_labels,
    style="top",
    plotter_kwargs=None,
    lineid_kwargs=None,
):
    """Use the lineid_plot package to plot user specified lines on a plotter's spectrum.

    Parameters
    ----------
    plotter : SDECPlotter
        Plotting interface for Spectral element DEComposition (SDEC) Plot.
    wavelengths : list of float
        wavelength values of the lines to be plotted
    line_labels : list of str
        the labels to be used for the lines at the wavelength positions
    style : {'top', 'inside', 'along spectrum'}, optional
        Preset styles for the lines markers. (default: 'top')
        If 'top', the lines will be plotted above the top axes.
        If 'inside', the lines will be plotted inside the axes.
        If 'along spectrum', the lines will be plotted along the spectrum.
    plotter_kwargs : dict, optional
        kwargs passed to the plotter's generate_plot_mpl method, by default None
    lineid_kwargs : dict, optional
        kwargs passed to the lineid_plot.plot_line_ids method, by default None

    Returns
    -------
    matplotlib.axes._subplots.AxesSubplot
        axes from generate_plot_mpl with the lines plotted.
    """

    if plotter_kwargs is None:
        plotter_kwargs = {}

    if lineid_kwargs is None:
        lineid_kwargs = {}

    ax = plotter.generate_plot_mpl(**plotter_kwargs)

    waves = plotter.plot_wavelength
    spec = plotter.modeled_spectrum_luminosity

    if style == "top" or style == None:
        fig, ax = lineid_plot.plot_line_ids(
            waves, spec, wavelengths, line_labels, ax=ax, **lineid_kwargs
        )

    elif style == "inside":
        box_loc = ax.transData.inverted().transform(
            ax.transAxes.transform((0, 0.9))
        )[1]
        arrow_loc = ax.transData.inverted().transform(
            ax.transAxes.transform((0, 0.8))
        )[1]

        fig, ax = lineid_plot.plot_line_ids(
            waves,
            spec,
            wavelengths,
            line_labels,
            ax=ax,
            arrow_tip=arrow_loc,
            box_loc=box_loc,
            **lineid_kwargs,
        )

    elif style == "along spectrum":
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
            wavelengths, waves.value[::-1], spec.value[::-1]
        )  # uses np.interp needs to monotonically inc

        arrow_tips = fluxes + marker_len
        box_locs = fluxes + marker_len + arrow_len

        fig, ax = lineid_plot.plot_line_ids(
            waves,
            spec,
            wavelengths,
            line_labels,
            ax=ax,
            arrow_tip=list(arrow_tips),
            box_loc=list(box_locs),
            **lineid_kwargs,
        )
    else:
        raise ValueError(
            "style must be one of 'top', 'inside', or 'along spectrum'"
        )

    return fig, ax
