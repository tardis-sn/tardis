import numpy as np
import numexpr as ne
from astropy import constants as const

class BlackBodySimpleSource(object):
    """
    Simple packet source that generates packets for the Montecarlo part.
    This uses the algorithm described in  Bjorkman & Wood 2001 (page 4) which
    references Carter & Cashwell 1975:

    First, generate a uniform random number, :math:`\\xi_0`, in
    the range 0 to 1, and determine the minimum value of $l$,
    $l_{\\rm min}$, that satisfies the condition
    %
    \\begin{equation}
    \\sum_{i=1}^{l} i^{-4} \\ge {{\\pi^4}\\over{90}} \\m_0 \\; .
    \\end{equation}
    %
    Next obtain four additional uniform random numbers (in the range 0 to 1),
    $\\xi_1$, $\\xi_2$, $\\xi_3$, and $\\xi_4$.  Finally, the packet
    frequency is given by
    %
    \\begin{equation}
    x = -\\ln{(\\xi_1\\xi_2\\xi_3\\xi_4)}/l_{\\rm min}\\; .
    \\end{equation}

    where :math:`x=h\\nu/kT`

    """
    def __init__(self, seed, l_samples=100000):
        np.random.seed=seed
        self.l_samples = l_samples
        self.l_array = np.cumsum(np.arange(1, l_samples)**-4)
        self.l_coef = np.pi**4 / 90.0

    def create_packets(self, T, no_of_packets):

        xis = np.random.random((5, no_of_packets))
        l = self.l_array.searchsorted(xis[0]*self.l_coef) + 1.
        if np.any(l > self.l_samples-10): 1/0
        xis_prod = np.prod(xis[1:], 0)
        x = ne.evaluate('-log(xis_prod)/l')

        return x * (const.k_B.cgs.value * T) / const.h.cgs.value
