Packet Sampling
===============


In the simulation, we use spherical symmetry where the ejecta is divided into multiple spherical shells(see :doc:`../setup/model`). To find the mass of each shell we multiply the mass fractions of each of the elements and isotopes with the density. We can set the timespace for the simulation
to be between t\ :sub:`start`\  and t\ :sub:`end`\  using either a linear or a logarithmic scale.


We use the `radioactivedecay <https://radioactivedecay.github.io>` python package to calculate the number of decays in each channel for each shell and at every timestep for each isotope.
The composition is also updated after each timestep following :cite:`Guttman2024`.
Then to calculate the decay energy of each channel, we multiply the total number of decays with the energy of each channel.
We give weight to the decay energy of each channel of each isotope to sample the packets which takes into account the spatial distribution
of radiocative isotopes present in different shells and the evolution of the composition over time.
