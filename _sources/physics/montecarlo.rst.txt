.. _montecarlo:

*********************
Radiative Monte Carlo
*********************

.. :currentmodule:: tardis.montecarlo_multizone

We assume the supernova is in homologous expansion and we assume that
there is photon injection from an inner boundary and we assume an outer boundary
where these photons can leave the simulation again.


Montecarlo packets
^^^^^^^^^^^^^^^^^^

Montecarlo packets have a frequency, angle (:math:`\mu=\cos{\theta}` )
and an energy: :math:`P(\nu, \mu, \epsilon)`. A large number :math:`n` are generated
at the start of each Montecarlo run in TARDIS by drawing :math:`\nu` from a
black-body distribution distribution with a :math:`\mu` being drawn from a
:math:`\sqrt[2]{z}` distribution and an
energy that is :math:`1/n`.

These packets are then launched into the simulation and there are two possible
outcomes for each packet.

 * **packet leaves the simulation through the outer bound** and count towards the
    spectrum
 * **packet leaves the simulation through the inner boundary** and are lost
    (reabsorbed)

The packets can gain and loose energy throughout the simulation. If these packets
scatter (through the various mechanisms), we transform into the comoving frame at
that position :math:`v(r) = r / t_\textrm{explosion}`, where :math:`t_explosion`
is the time since explosion with the doppler factor
:math:`d_\textrm{rest \rightarrow comoving}` then change direction from
:math:`\mu_\textrm{in}` to :math:`mu_\textrm{out}` and transform back and
change the energy accordingly:

.. math::
    d_\textrm{rest \rightarrow comoving} = 1 - \mu_\textrm{in} v(r) / c \\
    d_\textrm{rest \rightarrow comoving} = 1 - \mu_\textrm{out} v(r) / c \\
    E_\textrm{after scatter} = E_\textrm{before scatter} \times \frac{1 - \mu_\textrm{in} v(r) / c}{1 - \mu_\textrm{out} v(r) / c}

This way the montecarlo packets can gain or loose energy in the simulation:

.. plot:: physics/pyplot/plot_mu_in_out_packet.py




The spectrum is then generated as a weighted histogram. For each bin with edges
:math:`\nu_{n}, \nu_{n+1}` we get all the packets that left the simulation through
the outer boundary and add up their remaining energies.



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

Algorithm Flowchart
^^^^^^^^^^^^^^^^^^^

.. graphviz::

  digraph g{
    a -> b -> c
    c -> d [label="d_inner or \nd_outer"]
    c -> e [label="d_line"]
    d -> f [label="yes"]
    d -> g [label="no"]
    g -> a
    e -> a [label="no"]
    e -> h [label="yes"]
    h -> a
    a [label="We have a packet.",shape=box,fillcolor="white",style="filled,rounded"];
    b [label="Calculate\nd_line, d_electron,\nd_inner and d_outer.",shape=box,fillcolor="white",style="filled,rounded"];
    c [label="Which distance\nis smallest?", shape="diamond", fillcolor="white", style="filled"]
    d [label="Are we leaving\nsimulation area?", shape="diamond", fillcolor="white", style="filled"]
    e [label="Does the\npacket interact?", shape="diamond", fillcolor="white", style="filled"]
    f [label="Packet is re-absorbed\nor emitted.\nThis ends the loop.", shape="box", fillcolor="white", style="filled,rounded"]
    g [label="Update line\nprobabilities.", shape="box", fillcolor="white", style="filled,rounded"]
    h [label="New random direction,\nupdated energy,\nmoving packet to current position,\nupdating event random number.", shape="box", fillcolor="white", style="filled,rounded"]
  }
