******************************
Modeling Supernovae Explosions
******************************


.. _expansion:

Homologous Expansion
====================

TARDIS models supernovae as expanding homologously. This means that at the beginning of the explosion, the supernova starts at a single point and proceeds to expand radially outward such that the ratio of the velocity of the ejecta to the distance from the ejecta to the supernova's center is uniform throughout the supernova. As an example, if the outer edge of the ejecta moves outward at some velocity :math:`v_{outer}`, the velocity of the ejecta half way between the outer edge and the center would be :math:`v_{outer}/2`. The animation below demonstrates this type of expansion.

TARDIS simulates radiative transfer between an inner boudary (inside of which the supernova is modeled as a blackbody-- the so-called "photosphere") and the outer boundary of the supernova (as explained in :ref:`propogation`). The velocity of the inner boundary :math:`v_{inner}` and the velocity of the outer boundary :math:`v_{outer}` is supplied in the configuration file, as well as the time after the explosion for which TARDIS is calculating the spectrum (:math:`t_{explosion}`). The radii of the inner and outer boundaries are therefore calcuated by :math:`r_{inner}=v_{inner}*t_{explosion}` and :math:`r_{outer}=v_{outer}*t_{explosion}`. Plasma a distance :math:`r` from the center of the supernova would then be traveling outward at a speed :math:`v_{plasma}=\frac{r}{r_{outer}}v_{outer}`. This is also shown in the animation.

Additionally, TARDIS divides the space between the inner and outer computational boundaries into cells-- radial shells for which the plasma state is (spacially) constant. In the animation, 6 cells are shown, being divided by the light blue lines. As TARDIS is a time-independent code which calculates the spectra at an instant in time, the radii boundaries (either of the computational domain or of the cells) do not chage throughout the simulation.


**ANIMATION**


.. _referenceframes:

Reference Frames
================

In TARDIS, two reference frames are of particular importance: the lab frame and the co-moving frame. In the lab frame, the center of the supernova is at rest-- for example, the animation above is shown in the lab frame. This is the frame for which the spectra are calculated.

The co-moving frame at some point in the supernova, however, has the plasma at that point be at rest. This is the frame of reference "according to the plasma."

If a photon is released from the photosphere with a frequency :math:`\nu_{lab}` in the lab frame, the doppler effect says that, in the co-moving frame at a distance :math:`r` from the center of the supernova, the photon's frequency is shifted to

.. math::
    \nu_{co-moving} = stuff.
    
where :math:`\beta = stuff` and :math:`v_{plasma}` is the radial velocity of the plasma.

**What if the photon is moving the other direction? Mention above about the photon moving outwards**





NOTES:

Do something like
-----------------

Numerical and physical events (mention TARDIS calculates distances)
        Distance to next cell
        Physical interations
            Talk about optical depth
            And the accumulation of optical depth
            And what happens after an interaction
        Example situations