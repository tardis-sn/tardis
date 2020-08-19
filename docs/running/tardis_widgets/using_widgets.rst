********************
Using TARDIS Widgets
********************

This page describes what each TARDIS Widget has to offer and how you can make
the best use of it. If you're looking for the code to generate widgets, head
over to `Generating TARDIS Widgets <generating_widgets>`_ section and see the
notebook in action.

Currently, TARDIS supports the following widgets:

Shell Info Widget
#################

This widget allows you to explore the chemical abundances in different shells
of a Supernova ejecta.

.. image:: images/shell-info-widget-demo.gif
    :alt: Demo of Shell Info Widget

It consists of four interlinked tables - clicking on any row in a table,
populates data in the table(s) right to it. All these 4 tables, in order, has
following information:

1. **Shells Data** - Radiative temperature and Dilution Factor (W) of each shell
(computational partitions) of a Supernova ejecta. Shell numbers in ascending 
order represent inner to outer i.e. innermost shell is numbered as 1 and 
outermost shell is numbered as the last number (say 20).

2. **Element Abundances** - Fractional abundance of each element present in the
selected shell.

3. **Ion Abundances** - Fractional abundance of each ion (species) of the
selected element present in the selected shell. 

4. **Level Abundances** - Fractional abundance of each level of the selected
ion and element in the selected shell.

.. any more info about ion no.s (0 means no charge and hence I) and level no.s

Line Info Widget
################

This widget lets you explore line interactions present in the spectrum of a
simulation model.

.. image:: images/shell-info-widget-demo.gif
    :alt: Demo of Line Info Widget

You can select any wavelength range in the spectrum interactively to display a 
table of species abundances (fraction of packets interacting). Using 
toggle buttons, you can specify whether to filter the selected range by
emitted or absorbed wavelengths of packets. Clicking on a row in species
abundances table further reveals packet counts for each last line 
interaction which can be grouped by excitation lines, de-excitation lines
or both, using the dropdown.

.. both tables and group by kinda self-explanatory, if required, can explain
what filter by actually does

Interacting with Spectrum
=========================
The spectrum in line info widget is an interactive figure made using plotly,
there are several things you can do with it:

Making Selection
----------------
By default, the box selection is enabled so you just need to click and drag on
figure and the pink colored selection box will appear. If you're wondering
why it has full height, it's because we're only interested in controlling
widths to specify range on X-axis i.e. wavelength range. 

.. image:: images/spectrum_selection.gif

After making a selection, if you need to resize the selection box (say make it
narrower), simply redraw a new selection box over the older one.


Using Rangesilder
-----------------
Since spectra are usually very long, so it makes sense to focus (horizontally 
zoom) only over a particular range on X-axis - this is what rangeslider (the 
long bar just below the figure) allows you to do!

.. image:: images/spectrum_rangeslider.gif

Either you can **slide** the focused range by clicking and dragging it or you 
can **resize** it by dragging the handles (vertical bars) at its edges.

Using other options in Modebar
------------------------------
If you take your mouse to the top right corner of figure, you will see a Modebar
with multiple options. The default option when line info widget first displays
is ***Box Select*** - the dotted square icon. You can click on other options
like ***Zoom** (magnifying glass icon), to do a rectangle zoom which maybe
helpful to focus on a feature in spectrum. You can always revert back to 

.. image:: images/spectrum_modebar.gif

There are also several other options in the modebar which we have not explained
(because they are not very relevant) but you're free to explore them as long as
you remember to click back on the Box Select option for making selections on
spectrum.

.. Toggle legend
