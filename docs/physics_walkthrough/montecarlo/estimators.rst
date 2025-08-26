.. _estimators:

***********************
Volume-Based Estimators
***********************

Besides from just tracking the propagation of our packets, TARDIS also uses the Monte Carlo iteration is to determine useful information about the light traveling through the supernova (also called the radiation field). This information will eventually be used to help :doc:`update the plasma state <../update_and_conv/update_and_conv>` as well as :doc:`generate different types of spectra <../spectrum/index>`. We determine this information through volume-based estimators. The concept was originally developed by :cite:`Lucy1999` and successively refined by :cite:`Lucy1999a`, :cite:`Lucy2002` and :cite:`Lucy2003`.


Theory
======

:math:`J` and :math:`\bar \nu` Estimators
-----------------------------------------

Ordinarily, TARDIS is not concerned about the physical amount of time a packet spends traveling through the ejecta. Instead, we consider the "time of simulation" :math:`\Delta t` which is chosen to be the amount of time in which the photosphere emits the ensemble of packets (see :doc:`Energy Packet Initialization <initialization>`). When looking at the estimators, a slightly different interpretation of the packets is necessary. Here, we view the packets as not carrying a discrete amount of energy :math:`\varepsilon` that is emitted in a time interval :math:`\Delta t`, but as being a flow of energy that carries an energy :math:`\varepsilon` over a time :math:`\Delta t` -- that is, each packet is carrying a luminosity (energy per unit time) of :math:`L = \frac{\varepsilon}{\Delta t}`. Now, we can say that if a packet spends a time :math:`\delta t` in the supernova's ejecta, it contributes an energy of :math:`L\delta t= \frac{\varepsilon}{\Delta t}\delta t` into the radiation energy of the ejecta.

To account for the effects of the Monte Carlo packets on the ejecta, TARDIS uses the packets to first determine the average radiation energy density :math:`E` throughout each shell, where the energy density is the total radiation energy in the shell divided by the volume of the shell :math:`V=\frac{4}{3}\pi (r_\mathrm{outer}^3-r_\mathrm{inner}^3)`. Therefore, we add up the amount of energy each packet contributes to the radiation energy in that shell, and divide by the total volume of the shell:

.. math:: E=\frac{1}{V}\sum_i L_i\delta t_i=\frac{1}{V}\sum_i \frac{\varepsilon_i}{\Delta t}\delta t_i = \frac{1}{V\Delta t}\sum_i \varepsilon_i\delta t_i

where we sum over every Monte Carlo packet in the shell. Note that we are interested in the energy density in the co-moving frame (i.e. the energy density "according to the plasma," see :ref:`referenceframes`). Now, we note that the amount of time the Monte Carlo packet spends in a shell is :math:`\delta t = \frac{l_i}{c}` where :math:`l` is the distance that the packet travels through the shell. Thus, our estimator is

.. math:: E=\frac{1}{V\Delta t}\sum_i \varepsilon_i\frac{l_i}{c} = \frac{1}{cV\Delta t}\sum_i \varepsilon_i l_i.

Using this energy density, we can then calculate the mean radiation intensity :math:`J` in that shell using the relation :math:`J=\frac{c}{4\pi} E`, which gives us

.. math:: J=\frac{1}{4\pi V\Delta t}\sum_i \varepsilon_i l_i.

Since along any path the co-moving energy of the packet is continuously doppler shifted, we approximate this estimator using the co-moving energy at the beginning of the packet's path (theoretically, the midpoint of the path would be a better option. However, we use the beginning of the path for computational ease at a very small cost to the estimator's accuracy).

Next, we calculate the mean radiation frequency in each shell. For this, in each shell we add up the frequency of each packet weighted by the intensity they contribute to the shell (and then divide by the total intensity, as will be discussed below). Remembering that intensity is :math:`\frac{c}{4\pi}` times the energy density, and as before each packet contributes an energy of :math:`\frac{\varepsilon_i l_i}{c\Delta t}` and thus energy density of :math:`\frac{\varepsilon_i l_i}{cV\Delta t}` to its shell, we have that each packet contributes an intensity of :math:`\frac{\varepsilon_i l_i}{4\pi V\Delta t}` to its shell. So,

.. math:: \bar \nu = \sum_i \frac{\varepsilon_i l_i}{4\pi V \Delta t}  \nu_i = \frac{1}{4\pi V \Delta t}\sum_i \varepsilon_i \nu_i l_i

where once again the co-moving energy and frequency of each packet are taken at the beginning of the packet's path.

It is important to note that, confusingly, :math:`\bar \nu` is not truly the mean frequency (as can be seen by its units -- it has dimensions of intensity times frequency). Indeed, the true mean frequency is actually :math:`\frac{\bar \nu}{J}`.

.. note:: Both estimators take on a different value in each shell.

These estimators will be used in the :doc:`../update_and_conv/update_and_conv` step to help update the plasma state between iterations.


.. _j-blue-estimator:

:math:`J^b_{lu}` Estimator
--------------------------

Another estimator, called the ``j_blue`` or :math:`J^b_{lu}` estimator, is unlike the two previous estimators discussed. Instead of storing the mean intensity over the entire spectrum, it stores the intensity at a specific frequency. More precisely, since frequency is a continuum, it stores the intensity per unit frequency. In each shell, we record the intensity per unit frequency at the blue end (higher frequency end; this is where the ":math:`b`" superscript in :math:`J^b_{lu}` comes from) of each line transition -- that is, if a line transition :math:`l\rightarrow u` (from the lower energy level :math:`l` to the upper energy level :math:`u`, hence the :math:`lu` in :math:`J^b_{lu}`) has a frequency :math:`\nu_{lu}`, the mean intensity between :math:`\nu_{lu}` and :math:`\nu_{lu}+d\nu` is :math:`J^b_{lu}d\nu`. **This means that the** :math:`J^b_{lu}` **estimator contains values for each atomic line in each shell**, and is hence called a *line estimator*. Now, using our previous :math:`J` estimator, we have

.. math:: J^b_{lu}d\nu = \frac{1}{4\pi V\Delta t}\sum_i \varepsilon_i dl_i

where :math:`dl_i` is the infinitesimal distance that the packet travels while it has a co-moving frequency between :math:`\nu_{lu}` and :math:`\nu_{lu}+d\nu` (here, therefore, we are summing over all packets in a shell that achieve a co-moving frequency of :math:`\nu_{lu}` and thus come into resonance with the line transition :math:`l\rightarrow u` within that shell).

Now, say the packet with lab frequency :math:`\nu_\mathrm{lab}` has a co-moving frequency of :math:`\nu_{lu}` at a radius :math:`r_1` and propagation direction :math:`\mu_1`, and it has a co-moving frequency of :math:`\nu_{lu}+d\nu` at a radius :math:`r_2` and propagation direction :math:`\mu_2`. Then (see :ref:`referenceframes`):

.. math:: \nu_{lu}=\left(1-\frac{r_1\mu_1}{ct_\mathrm{explosion}}\right)\nu_\mathrm{lab}

and

.. math:: \nu_{lu}+d\nu=\left(1-\frac{r_2\mu_2}{ct_\mathrm{explosion}}\right)\nu_\mathrm{lab}.

But then subtracting, we get

.. math:: d\nu = (r_2\mu_2-r_1\mu_1)\frac{\nu_\mathrm{lab}}{ct_\mathrm{explosion}}=dl*\frac{\nu_\mathrm{lab}}{ct_\mathrm{explosion}}

(for the last equality, see :ref:`spherical-domain`).

But now inputting this into the equation for :math:`J^b_{lu}` (using :math:`\frac{dl_i}{d\nu}=\frac{ct_\mathrm{explosion}}{\nu_\mathrm{lab,i}}`), we get

.. math:: J^b_{lu} = \frac{ct_\mathrm{explosion}}{4\pi V\Delta t}\sum_i \frac{\varepsilon_i}{\nu_\mathrm{lab,i}}.


.. _edotlu:

:math:`\dot{E}_{lu}` Estimator
------------------------------

The final estimator we consider, like the ``j_blue`` estimator, is a line estimator, meaning it has contains values for each atomic line in each shell. It calculates the rate at which energy density interacts with a line transition :math:`l\rightarrow u`. The first step is to calculate the rate at which energy density resonates with some line in some shell. Each packet that comes into resonance with the transition :math:`l\rightarrow u` in a shell of volume :math:`V` contributes an energy density to that shell of :math:`\frac{\varepsilon}{V}` over a time :math:`\Delta t`, meaning the rate at which energy density resonates with the line is :math:`\sum_i \frac{\varepsilon_i}{\Delta t V} = \frac{1}{\Delta t V} \sum \varepsilon` where we are summing over all packets which come into resonance with the atomic line in some shell (as usual, this sum is done using the energies in the co-moving frame). Finally, this light then has a :math:`\left( 1- e^{-\tau_{lu}}\right)` probability of interacting with the line (where :math:`\tau_{lu}` is the Sobolev optical depth for the transition :math:`l\rightarrow u`, see :ref:`physical-interactions`), so the rate at which energy density is absorbed into the transition :math:`l\rightarrow u`, called the ``Edotlu`` estimator, is

.. math:: \dot{E}_{lu} = \frac{1}{\Delta t V} \left( 1- e^{-\tau_{lu}}\right) \sum_i \varepsilon_i.

Note that while one may expect us to merely add up the contributions of each packet that *interacts* with the line, eliminating the need for the :math:`\left( 1- e^{-\tau_{lu}}\right)` term, the former determines the desired quantity with more accuracy and less noise than the latter, (this is because it does not depend on the limited number of packets TARDIS uses, rather a theoretical equation, to determine how much of the light that comes into resonance with a line actually interacts with that line).


Implementation
==============

As previously discussed, a major component of each Monte Carlo iteration is the packet propagation process. During the packet propagation process this step, the :math:`J` and :math:`\bar \nu` estimators are updates every time a packet is moved to the next event location. Specifically, every time a packet is moved, :math:`\varepsilon l` is added to the "running total" :math:`J` estimator in the shell where the packet is, and :math:`\varepsilon \nu l` is added to the "running total" :math:`\bar\nu` estimator in the shell where the packet is (where :math:`l` is the distance the packet is moved, and :math:`\varepsilon` and :math:`\nu` are respectively the packet's co-moving energy and frequency at the beginning of the packet's path). The factor of :math:`\frac{1}{4\pi V\Delta t}`, for computational ease, is not attached to the estimators but is included during any calculations using these estimators, see :doc:`../update_and_conv/update_and_conv`.

Additionally, during the propagation process, every time a packet passes through a Sobolev point, meaning it reaches a co-moving frequency of :math:`\nu_{lu}` for some transition :math:`l\rightarrow u` and thus comes in resonance with an atomic line, the :math:`J^b_{lu}` for that atomic transition in the shell it is in is incremented by :math:`\frac{\varepsilon}{\nu_\mathrm{lab}}`, where :math:`\varepsilon` is the packet's energy and :math:`\nu_\mathrm{lab}` is the packet's lab-frame frequency. As before, for computational ease, the factor :math:`\frac{ct_\mathrm{explosion}}{4\pi V \Delta t}` is included in calculations involving the estimator (such as when `updating <../update_and_conv/update_and_conv.ipynb#updating-other-quantities>`_ :math:`J^b_{lu}` in the plasma or in the :ref:`formal integral <formal_integral>`). Similarly, when a packet passes through a Sobolev point, the :math:`\dot{E}_{lu}` for that atomic transition in the shell it is in is incremented by :math:`\varepsilon`, and once again, for computational ease, the term :math:`\frac{1}{\Delta t V} \left( 1- e^{-\tau_{lu}}\right)` is not included until calculations involving the estimator are performed (specifically in the :ref:`formal integral <formal_integral>`).

.. note:: Since the ``j_blue`` and ``Edotlu`` estimators are updated every time a packet comes into resonance with an atomic line (not necessarily going through a line interaction), these estimators are equal to zero in some shell for a specific line if (and only if) no packet comes into resonance with that line within the shell.
