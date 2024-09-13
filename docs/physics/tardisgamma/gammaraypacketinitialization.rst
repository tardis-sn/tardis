*********************
Tardis :math:`\gamma`
*********************

Packet Initialization
=====================

The gamma ray portion of tardis also uses packets of photons, which are given the following properties:
- **Packet Energy:** This is the energy of the packet in the comoving frame and is equal to the total energy divided by the number of packets.
- **Packet Frequency:** This is the frequency of the packet in the comoving frame and is equal to the packet energy divided by Planck's constant, h.
- **Packet Direction:** This is the direction of the packet and is made up of two angles, :math:`\theta` and :math:`\phi`. 
:math:`\theta` is a polar angle between 0 and :math:`\pi` and :math:`\phi` is an azimuth angle between 0 and 2:math:`\pi`.
We randomly sample these angles by using a random numbers z\ :sub:`1`| and z\ :sub:`2`\.
We then use the equations :math:`\cos{\theta}`= 1-2z\ :sub:`1`\ and :math:`\phi` = 2:math:`\pi`z\ :sub:`2`\ to find the angles.
- **Packet Times** This is the time that the packet starts propagating the ejecta.
- **Packet Location** This is the starting point of the packet within the spherical shells of the ejecta.
This location is given by the equation, v = [zv\ :sub:`inner`\\ :sup:`3`\ + (1-z)v\ :sub:`outer`\\ :sup:`3`\]\ :sup:`1/3`\
where v\ :sub:`inner`\ and v\ :sub:`outer`\ are the inner and outer velocities of the shell and z is a random number between [0,1).
Then to get the radial position, r, we multiply this velocity by the packet time. To get the Cartesian coordinates we use the equations:
x = r:math:`\sin{\theta}\cos{\phi}`  y = r:math:`\sin{\theta}\cos{\phi}`  z = r:math:`\cos{\theta}`
