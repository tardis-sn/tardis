.. _two_photon_archive:

Archived two-photon continuum calculations
------------------------------------------

This page preserves the two-photon-decay calculations that were removed with
the legacy continuum plasma properties. They are documentation only: the
classic plasma workflow no longer assembles or evaluates these properties.

Inputs and transition weight
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each tabulated two-photon decay, the legacy implementation used the
atomic-data fields :math:`A_{ul}`, :math:`\nu_0`, and the fit coefficients
:math:`\alpha`, :math:`\beta`, and :math:`\gamma`. Here :math:`\nu_0` is the
frequency of the corresponding normal line transition. The unnormalised
macro-atom deactivation weight was

.. math::

    p_{2\gamma} = A_{ul} h \nu_0.

The same scalar weight was copied to every model shell. It was registered as a
deactivation channel labelled ``two-photon``; the legacy code did not model
internal two-photon transitions.

The corresponding implementation was:

.. code-block:: python

    no_shells = len(density)
    p_two_phot = two_photon_data.A_ul * two_photon_data.nu0 * H
    p_two_phot = pd.concat([p_two_phot] * no_shells, axis=1)

``H`` was the Planck constant in CGS units. The remaining legacy code attached
the macro-atom source and destination indices and labelled the destination
channel ``two-photon``.

Spectral shape
~~~~~~~~~~~~~~

For a photon frequency :math:`\nu \in [0, \nu_0]`, the reduced frequency is

.. math::

    y = \frac{\nu}{\nu_0}.

The legacy implementation evaluated the fitted spectral shape described by
Nussbaumer & Schmutz (1984):

.. math::

    A(y) = y(1-y)\left[1 - \left(4y(1-y)\right)^\gamma\right]
           + \alpha\left[y(1-y)\right]^\beta
             \left(4y(1-y)\right)^\gamma,

    j_\nu = y A(y).

:math:`j_\nu` was intentionally left unnormalised because it was used only to
form relative emission probabilities.

The implementation of the fitted emissivity was:

.. code-block:: python

    def calculate_j_nu(y, alpha, beta, gamma):
        ay = y * (1 - y) * (1 - (4 * y * (1 - y)) ** gamma)
        ay += alpha * (y * (1 - y)) ** beta * (4 * y * (1 - y)) ** gamma
        return ay * y

Tabulated CDF and sampling
~~~~~~~~~~~~~~~~~~~~~~~~~~

For every transition, the implementation tabulated 500 uniformly spaced
frequencies from :math:`0` through :math:`\nu_0`. It constructed the
cumulative distribution by trapezoidal integration,

.. math::

    C(\nu_i) = \frac{\displaystyle\int_0^{\nu_i} j_\nu\,d\nu}
                     {\displaystyle\int_0^{\nu_0} j_\nu\,d\nu}.

To sample an emitted frequency, it drew :math:`z \sim U(0, 1)`, selected the
first tabulated CDF value greater than :math:`z`, and linearly interpolated
between that point and the preceding point. Given adjacent samples
:math:`(\nu_{i-1}, C_{i-1})` and :math:`(\nu_i, C_i)`, this was

.. math::

    \nu = \nu_i - \frac{C_i-z}{C_i-C_{i-1}}
                    \left(\nu_i-\nu_{i-1}\right).

The CDF construction and sampler were implemented as follows:

.. code-block:: python

    bins = 500
    nu = np.linspace(0.0, row.nu0, bins)
    y = nu / row.nu0
    j_nu = calculate_j_nu(y, row.alpha, row.beta, row.gamma)

    cdf = np.zeros_like(nu)
    cdf[1:] = numba_cumulative_trapezoid(j_nu, nu)
    cdf /= cdf[-1]

    zrand = np.random.random()
    idx = np.searchsorted(cdf, zrand, side="right")
    nu_sample = nu[idx] - (cdf[idx] - zrand) / (
        cdf[idx] - cdf[idx - 1]
    ) * (nu[idx] - nu[idx - 1])

The result is a sample from the tabulated two-photon emission spectrum. This
algorithm is retained here for scientific traceability and is not an available
classic-plasma API.
