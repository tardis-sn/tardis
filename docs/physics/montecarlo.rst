.. _montecarlo:

*********************
Radiative Monte Carlo
*********************

.. :currentmodule:: tardis.montecarlo_multizone


Montecarlo Geometry
^^^^^^^^^^^^^^^^^^^

Before any packet action is performed we calculate four different distances
 ( :math:`d_\textrm{inner}, d_\textrm{outer}, d_\textrm{line}, d_{\textrm{e}^{-}}` )

The calculations for the distance to the outer boundary:

.. image:: ../graphics/d_outer.png
    :width: 400

The calculations for the distance to the inner boundary:

.. image:: ../graphics/d_inner.png
    :width: 400




Radiationfield estimators
^^^^^^^^^^^^^^^^^^^^^^^^^

During the monte-carlo run we collect two estimators for the radiation field:

.. math::

    J_\textrm{estimator} &= \sum{\epsilon l}\\
    \bar{\nu}_\textrm{estimator} &=  \sum{\epsilon \nu l},

where :math:`\epsilon, \nu` are comoving energy and comoving frequency of a packet respectively.

To calculate the temperature and dilution factor we first calculate the mean intensity in each cell
( :math:`J = \frac{1}{4\pi\, \Delta t\, V} J_\textrm{estimator}` )., :cite:`2003A&A...403..261L`.

The weighted mean frequency is used to obtain the radiation temperature. Specifically, the radiation temperature is chosen as the 
temperature of a black body that has the same weighted mean frequency as has been computed in the simulation. Accordingly,

.. math::

    \frac{h \bar{\nu}}{k_{B} T_{R}} = \frac{h}{k_{B} T_{R}} \frac{\bar{\nu}_\textrm{estimator}}{J_\textrm{estimator}} 
      = 24 \zeta(5) \frac{15}{\pi^4},

where the evaluation comes from the mean value of

.. math::

    \bar{x} = \frac{ \int_0^{\infty} x^4 / (\exp{x} - 1)dx}{\int_0^{\infty} x^3 / (\exp{x} - 1)dx} =
    24 \zeta(5) \frac{15}{\pi^4} = 3.8322\dots

and so

.. math::

    T_{R} &= \frac{1}{\bar{x}} \frac{h}{k_{B}} \frac{\bar{\nu}_\textrm{estimator}}{J_\textrm{estimator}} \\
    &= 0.260945 \frac{h}{k_{B}} \frac{\bar{\nu}_\textrm{estimator}}{J_\textrm{estimator}}.

With the radiation temperature known, we can then obtain our estimate for for the dilution factor. Our radiation field model in the 
nebular approximation is

.. math::

    J = W B(T_{R}) = W \frac{\sigma_{SB}}{\pi} T_{R}^4,

i.e. a dilute blackbody. Therefore we use our value of the mean intensity derrived from the estimator (above) to obtain the 
dilution factor

.. math::

    W = \frac{\pi J}{\sigma_{SB} T_{R}^4} = \frac{1}{4\sigma_{SB} T_{R}^4\, \Delta t\, V} J_\textrm{estimator}.

There endeth the lesson.
