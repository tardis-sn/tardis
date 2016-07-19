********************************
Direct source integration method
********************************

.. note::

  The following is provisional information, describing a code path that is right now active in the downbranch scheme.


One way to increase the speed of the monte carlo procedure is by improving final method of generating the spectra so that the quality of spectra produced by a given amount of packets is increased. This is the goal of the source integration scheme of :cite `Lucy99b`, which replaces the simple binning of the escaping packets with a method based on the formal integral of the emergent intensity.

The procedure starts with a monte carlo line absorption rate estimator:

.. math::

   \dot E_{lu} = \frac{1}{\Delta t V} \left( 1- e^{-\tau_lu}\right) \sum \epsilon

where the sum is over all the packages in a given shell that come into resonance with the transition :math:`u \rightarrow l` during the monte carlo run, :math:`\epsilon` is the energy of one such packet, and :math:`\tau_{lu}` the optical depth of the line. The sum of estimator is implemented in the c code as `increment_Edotlu_estimator` and the prefactor is calculated in the `postprocess` function called at the end of `montecarlo_radial1d` in `montecarlo.pyx` if the `last_run` argument is set to `True`. Right now indicating the last run is done in `legacy_run_simulation` by way of a hard coded value. 

After the final monte carlo step, a level absorption estimator is calculated, given by:

.. math::

   \dot E_u = \sum_{i < u} \dot E_{lu}

that is, by summing all the line absorption estimators below the curently selected level. In the code this is done with a bit of pandas magic in `make_source_function` found in `simulation/base.py`. By creating a dataframe using the estimated :math:`\dot E_u` and appyling to it a copy of the index of the atomic data, the sum above can be done with to a `groupby` operation following by a sum of the result. 

The source function for each line can then be derived from the relation

.. math::

   \left( 1- e^{-\tau_lu}\right) S_{ul} = \frac{\lambda_{ul} t}{4 \pi} q_{ul} \dot E_u

where :math:`\lambda_{ul}` is the wavelength of each line  :math:`u \rightarrow l`, and :math:`q_{ul}` is the corresponding branching ratio. The attenuating factor is kept on the left hand side because it is the product of the two that will appear in later formulae. The product on the right hand side is also evaluated in `make_source_function`. 

Having thus produced attenuated source functions from our Monte Carlo run, we move on to using this to calculate the emerging intensity and finally the luminosity per wavelength. The final integral is given as 

.. math::

   L_\nu  = 8 \pi^2 \int_0^\infty I_\nu (p) p dp

where :math:`p` is the impact parameter of a ray trough the supernova envelope that reaches the distant observer, and :math:`I_\nu (p)` is the intensity along one such ray, given by recursing trough the list of attenuated source functions from the blue to the red and adding up contributions. The relation linking the intensity before the k:th transition :math:`u \rightarrow l` to the intensity after is 

.. math::

   I_k^r = I_k^b e^{-\tau_k} + \left( 1- e^{-\tau_k}\right) S_{k}

where the superscripts are crucial, with :math:`r` and :math:`b` referencing the red and blue sides of the k:th transition respectively. To go from the red side of a line to the blue side of the next we can either ignore continuum sources of opacity, in which case

.. math:: 

   I_{k+1}^b = I_k^r

or include them, then requiring we perform

.. math:: 

   I_{k+1}^b = I_k^r + \Delta \tau_k \left[ \frac 1 2(J_k^r + J_{k+1}^b) - I_k^r  \right]

The starting condition for the blue to red side transition is either :math:`I_0^r = 0` for the case that the impact parameter is greater than the radius if the photosphere, or :math:`I_0^r = B_\nu(T)` if the impact parameter is less than the radius of the photosphere. 

.. note::

   Currently the code does not perform the steps necessary to include continuum sources of opacity.

We seek to integrate all emissions at a certain wavelength :math:`\nu` along a ray at some specific impact parameter :math:`p'`. Because the SNE is expanding homologously, along any ray parallel to the line of sight the Doppler effect will sift a range of co-moving frequencies :math:`\nu_0` into the desired observer frame frequency :math:`\nu`.

To find which lines to include when recursing on the line list we need to find the maximum Doppler shift along a given ray. At any point, the Doppler effect in our coordinates is

.. math::

   \nu = \nu_0 \left( 1 + \beta \mu \right)

where :math:`\beta = \frac v c`, and :math:`\mu = \cos \theta`. Here :math:`\theta` is the angle between the radial direction and the ray to the observer, as shown in the figure below. Because we are in the homologous expansion regime :math:`v = \frac r t`. Solving for :math:`\nu_0` in the above gives the relation we seek, but we require an expression for :math:`\mu`. Examining the figure, we see that for positive :math:`z` the angle :math:`\theta_2` can be related to the :math:`z` coordinate of the point C by

.. math::

   \cos \theta_2 = \frac{z_c}{r} = \mu 

.. image:: https://i.imgur.com/WwVHp5c.png

and in turn :math:`z_c` can be given as 

.. math::

   z_c = \sqrt{r_c^2 + p_c^2}

where the subscripts indicate the value at point C. By symmetry the intersection point for negative :math:`z` is simply :math:`-z_c`.

Using the expression for :math:`\mu`, :math:`\beta` above leads to the dependence on :math:`r` cancelling, and solving for :math:`\nu_0` gives

.. math::

   \nu_0 = \frac{\nu}{1 + \frac{z}{ct}}

For any given shell and impact parameter we can thus find the maximum and minimum co-moving frequency that will give the specified lab frame frequency if we know the intersection points of the ray with correct impact parameter, and this we find easily given the impact parameter and the radius of the shell.

The integrator function proceeds as follows: first set up a grid of relative impact parameters `ps`. Then take all the shell radii from largest to smallest and put them in relative units by dividing with the largest radius. Using these two we calculate the :math:`z` coordinate of the crossing points if the various impact parameters, directly yielding the positive :math:`z_c` in the upper triangular matrix `z_crossings`. I also include the normalization factor :math:`ct`.

Because the recursions have different starting conditions, we then split the crossings and impact parameters into an 'inner' and an 'outer' part, defined by whether the impact parameter a crossing corresponds to is greater or smaller than the innermost shell radius `R_min_rel`.

Then we simply iterate over all the frequencies we want, and for each frequency over both sections of impact parameters, in each section recursing over selections of the line lists derived by from the crossing points.
