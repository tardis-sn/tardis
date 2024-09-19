*********************
Tardis :math:`\gamma`
*********************

Packet Initialization
=====================

The gamma ray portion of tardis also uses the idea of packets of photons from :cite:`AbbottLucy85` and :cite:`Lucy99` (see :doc:`initialization`)
These packets are given an energy equal to the total energy divided by the number of packets and a frequency that is equal to the packet energy divided by Planck's constant, h.

The packets are also given a direction made up of two angles, :math:`\theta` and :math:`\phi`, where  :math:`\theta` is a polar angle between 0 and :math:`\pi` and :math:`\phi` is an azimuth angle between 0 and :math:`2\pi`.
We sample these angles using the following equations :cite:`CarterCashwell75`:

:math:`\cos{\theta}` = 1-2 z \ :sub:`1`\

:math:`\phi` = :math:`2\pi` z\ :sub:`2`\

The packets are also given a time and location where they start propagating the ejecta. We use the following equations to give the starting location:

v = [zv\ :sup:`3`\ :sub:`inner`\ + (1-z)v\ :sup:`3`\ :sub:`outer`\]\ :sup:`1/3`\

where v\ :sub:`inner`\  and v\ :sub:`outer`\  are the inner and outer velocities of the shell and z is a random number between [0,1).

Then to get the radial position, r, we multiply this velocity by the packet time.

Finally, to get the Cartesian coordinates we use the equations:

x = r :math:`\sin{\theta}\cos{\phi}`

y = r :math:`\sin{\theta}\cos{\phi}`

z = r :math:`\cos{\theta}`


