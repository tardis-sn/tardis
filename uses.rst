**************************
What can I do with TARDIS?
**************************

TARDIS is designed to carry out calculations of synthetic supernova spectra (see Kerzendorf & Sim 2014).

A TARDIS calculation requires

1) a model for the supernova ejecta (density and composition distribution);
2) that parameters determining the physical and numerical properties of the code are given.

These inputs are specified via the input file (yaml file). As a starting point, examine the input file in the example download (above). The example file sets up a very simple model in which the density distribution is based on the well-known W7 model (Nomoto et al. 1984; fit following Branch et al. 1985) and in which uniform abundances are specified. To get you started, here is a short list of some settings in that example that can be experimented with to get a feel for using TARDIS. Please contact us if you have questions or would like information on more advanced possibilities!

Change the abundances
=====================

In the example file, you can alter the elemental abundances (mass fractions; specified in the "abundances" section). Use this to explore the sensitivity of features to composition.

Change the luminosity
=====================

You can alter the (emergent) luminosity and see how this affects the synthetic spectrum. Do this by varying the "luminosity requested" in the "supernova" section.

Change the epoch
================

In the example file, you can change the epoch for which the synthetic spectrum is calculated (change the value of "time explosion"; specified in the "supernova" section). When doing this you might also change the inner boundary velocity ("start" value in the "velocity" section), and probably the luminosity (see above)!

Experiment with the treatment of line opacity
=============================================

In the "plasma" section you can change the "line interaction type" between "scatter", "downbranch" and "macroatom" - investigate how important fluorescence is in your synthetic spectrum. (See Kerzendorf & Sim and references therein for the meaning of these settings.)



