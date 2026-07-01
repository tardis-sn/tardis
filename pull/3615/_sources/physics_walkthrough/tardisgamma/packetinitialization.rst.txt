Packet Initialization
=====================

The gamma ray portion of TARDIS also uses the idea of packets of photons from :cite:`Abbott1985` and :cite:`Lucy1999` (see :doc:`../montecarlo/initialization`)
These packets are given an energy in the comoving frame which is equal to the total energy divided by the number of packets. They are also given a  frequency that is equal to the packet energy divided by Planck's constant, h.

The packets are also given a direction made up of two angles, :math:`\theta` and :math:`\phi`, where  :math:`\theta` is a polar angle between 0 and :math:`\pi` and :math:`\phi` is an azimuth angle between 0 and :math:`2\pi`.
We sample these angles using the following equations :cite:`Carter1975`:

.. math::

    \cos{\theta} = 1-2 z_1

    \phi = 2\pi z_2

The packets are also given a time and location where they start propagating the ejecta. We use the following equations to give the starting location:

.. math::
    v = \left[zv_{\text{inner}}^3 + (1-z)v_{\text{inner}}^3\right]^{1/3}

where v\ :sub:`inner`\  and v\ :sub:`outer`\  are the inner and outer velocities of the shell and z is a random number between [0,1).

Then to get the radial position, r, we multiply this velocity by the packet time.

Finally, to get the Cartesian coordinates we use the equations:

.. math::
    x = r\sin{\theta}\cos{\phi}

    y = r\sin{\theta}\cos{\phi}

    z = r\cos{\theta}