****************************************
Spectrum Generation with Virtual Packets
****************************************

The main purpose of TARDIS is the generation of synthetic spectra. Currently,
two methods are implemented to calculate the spectrum during the main Monte
Carlo calculation. One follows the obvious approach of recording the properties
of all escaping Monte Carlo packets and binning their contributions in
frequency (or wavelength) space. This real packet spectrum will naturally
suffer from Monte Carlo noise and if one tries to improve its signal-to-noise
ratio, one immediately encounters a fundamental property of Monte Carlo
approaches. Since Monte Carlo processes are typically governed by Poisson
statistics, the level of stochastic fluctuations decreases only as :math:`\propto
N^{-\frac{1}{2}}`, with :math:`N` denoting the number of Monte Carlo
experiments. Thus, to decrease the noise by an order of magnitude 100 times
more experiments have to be performed. In the case of the real packet spectrum
this translates into using 100 times more packets, which would increase the
runtime by about the same factor.

.. note::

    More details about Monte Carlo errors and noise behaviour may be found in
    the standard literature, for example in :cite:`Kalos2008`.

It is difficult to avoid this fundamental behaviour of Monte Carlo techniques.
However, sophisticated Monte Carlo techniques exist, which make better use of
the computational resources. One such approach, which achieves a better noise
behaviour for the same computational costs, is implemented in TARDIS. It relies
on the concept of so-called virtual packets and goes back to the works by
:cite:`Long2002` and :cite:`Sim2005`.

Virtual Packets
===============
