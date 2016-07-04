********************************
Direct source integration method
********************************

.. note::

    The following is provisional information, describing a code path that is right now active in the downbranch scheme.


One way to increase the speed of the monte carlo procedure is by improving final method of generating the spectra so that the quality of spectra produced by a given amount of packets is increased. This is the goal of the source integration scheme of :cite `Lucy99b`, which replaces the simple binning of the escaping packets with a method based on the formal integral of the emergent intensity.

The procedure starts with a monte carlo line absorption rate estimator:
.. math::

    \dot E_{lu} = \frac{1}{\Delta t V} \left( 1- e^{-\tau_lu}\right) \sum \epsilon

where the sum is over all the packages in a given shell that come into resonance with the transition :math:`u \right l` during the monte carlo run, :math:`\epsilon` is the energy of one such packet, and :math:`\tau_{lu}` the optical depth of the line. The sum of estimator is implemented in the c code as `increment_Edotlu_estimator` and the prefactor is calculated in the `postprocess` function called at the end of `montecarlo_radial1d` in `montecarlo.pyx` if the `last_run` argument is set to `True`. Right now indicating the last run is done in `legacy_run_simulation` by way of a hard coded value. 

After the final monte carlo step, a level absorption estimator is calculated, given by:

.. math::

    \dot E_u = \sum_{i < u} \dot E_{lu}

that is, by summing all the line absorption estimators below the curently selected level. In the code this is done with a bit of pandas magic in `make_source_function` found in `simulation/base.py`. By creating a dataframe using the estimated :math:`\dot E_u` and appyling to it a copy of the index of the atomic data, the sum above can be done with to a `groupby` operation following by a sum of the result. 

The source function for each line can then be derived from the relation

.. math::
    \left( 1- e^{-\tau_lu}\right) S_{ul} = \frac{\lambda_{ul} t}{4 \pi} q_{ul} \dot E_u

where math::`\lambda_{ul}` is the wavelength of each line  :math:`u \right l`, and math::`q_{ul}` is the corresponding branching ratio. The attenuating factor is kept on the left hand side because it is the product of the two that will appear in later formulae. The product on the right hand side is also evaluated in `make_source_function`. 

