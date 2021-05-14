*********************************************
Monte Carlo Discretization --- Energy Packets
*********************************************

While it is instructive to think about tracking the propagation history of
photons when illustrating the basic idea behind Monte Carlo radiative transfer
techniques, there are important numerical reasons for using a different
discretization scheme. Instead of thinking in the photon picture, it brings
significant advantages to follow the idea of :cite:`Abbott1985` and
:cite:`Lucy1999` and consider parcels of radiant energy as the fundamental
building blocks of the Monte Carlo calculation. These basic Monte Carlo quanta
are commonly referred to as "energy packets" or simply "packets".

During a Monte Carlo calculation, a large number of packets, all with a certain
energy :math:`\varepsilon`, are created. In addition, each packet is associated
with a frequency. These assignments are performed in a manner which ensures
that the ensemble of packets represents the spectral energy distribution of the
radiation field (see :ref:`Propagation <propagation>`).

During the simulation, the energy of the packet remains constant in the local
co-moving frame (see :ref:`Reference Frames <referenceframes>` for
details about the lab and co-moving frames). This naturally ensures energy
conservation and constitutes the main advantage of this discretization scheme.
There is one side effect of this so-called indivisible packet energy scheme
which often causes confusion: Even during radiation-matter interactions the
packet energy is conserved in the co-moving frame (see :doc:`Propagation
<propagation>`). However, the frequency associated with a packet may change
(e.g. during non-resonant line interactions). As a consequence, packets may
represent a varying number of real photons during their lifetime.

.. note::
    The indivisible energy packet scheme does not require that all packets have
    the same energy. This is just a convenient and simple choice adopted in
    TARDIS.

