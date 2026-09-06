***********
Positronium
***********

Positronium (Ps) is the pairing of an electron and positron into a bound state analogous to a hydrogen atom (see `Wikipedia <https://en.wikipedia.org/wiki/Positronium>`_). This occurs when the relative motion between the positron and electron is small :cite:`Jauch1976` (:math:`\beta \rightarrow 0`, :math:`\beta=v/c`).

The mean lifetime of Ps is very short, on the order of nanoseconds. It is most likely to produce 2 or 3 gamma-ray photons depending on if it is para-Ps (2 photons of 511 keV each) or ortho-Ps (3 photons). The Ps types are dependent on the spin state of the electron and so occur in a 1/4 (para-Ps, S = 0, Ms = 0) to 3/4 (ortho-Ps, S = 1, Ms = −1, 0, 1) ratio. The triplet emission produces a continuum of photon energies described in :cite:`Ore1949`.

In supernovae, Ps can form when :math:`\beta`-decay occurs in the radioactive material produced as part of the explosion. The amount of Ps formed compared to immediate 2 photon annihilation in supernovae is unknown. In the galaxy, the value is determined to be :math:`0.94\pm0.04` :cite:`Milne2004`. Gamma-ray transport codes typically assume either no Ps formation (and thus all :math:`\beta`-decay releases 2 511 keV photons) or 100% Ps formation (and thus photons are released in pairs or triplets in a 1/4 to 3/4 ratio).

According to :cite:`Jauch1976` the density-dependent reciprocal lifetime :math:`\frac{1}{\tau_2}` (called a "cross-section") of para-Ps is

.. math:: 
    P \equiv \frac{1}{\tau_2} = r_0^2 \pi \rho

where :math:`r_0` is the classical electron radius and :math:`\rho` is the electron number density. Note that the units here appear to be inverse length rather than the expected inverse time.

For ortho-Ps, the cross-section is 

.. math:: 
    P \equiv \frac{1}{\tau_3} = \frac{4}{3} (\pi^2-9) \alpha r_0^2  \rho

where :math:`\alpha` is the fine structure constant. This makes the ortho-Ps cross-section (and associated lifetime) roughly 300 times smaller than that of para-Ps.