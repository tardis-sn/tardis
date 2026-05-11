import numpy as np
from numba import njit

from tardis.transport.montecarlo import (
    njit_dict_no_parallel,
)


@njit(**njit_dict_no_parallel)
def depressed_quartic(A: float, B: float, C: float, D: float, E: float, threshold_fac=6.0e-7):

    a = -3.0*B**2.0/(8.0*A**2.0) + C/A
    b = B**3.0/(8.0*A**3.0) - B*C/(2.0*A**2.0) + D/A
    c = -3.0*B**4.0/(256.0*A**4.0) + C*B**2.0/(16.0*A**3.0) - B*D/(4.0*A**2.0) + E/A

    p = -a**2.0/12.0 - c
    q = -a**3.0/108.0 + a*c/3.0 - b**2.0/8.0


    # Certain solutions will have (physical) complex values of r
    # but the complex component can cancel out in y
    q2p3_term = q**2.0/(4.0) + p**3.0/27.0
    q2p3_complex = q2p3_term + 0j
    r = -q/2 + np.sqrt(q2p3_complex)
    u = r**(1./3.)

    # Handle case where u=0
    if u == 0.0:
        y = -5.0/6.0 * a - q**(1.0/3.0)
    else:
        y = -5.0/6.0 * a + u - p/(3.0*u)

    # TODO: potentially handle imaginary parts here, but they should be close to zero
    # It's possible there's a better way to threshold these to ensure that the cancellations
    # between terms are as perfect as they can be, and that may involve working with the
    # complex components of this y term. For now though, I'm discarding the imaginary part
    # and doing a slightly messier threshold with the "aybw" term below instead.
    #y_real = np.real(y)
    #y_imag = np.imag(y)

    y = np.real(y)
    w = np.sqrt(a + 2.0*y)

    # Handle floating point error to prevent small neg numbers in sqrt
    aybw_thresh = threshold_fac * max(np.abs(3.0*a), np.abs(2.0*y), np.abs(2.0*b/w))

    # connor-mcclellan: note - we lose a lot more precision here than I expect.
    # Edge cases require a threshold as large as 1e-7 to recover the same distance
    # to line as homologous expansion, but I would have expected ~1e-15 would work
    # Somewhat clumsy syntax here is for numba compatibility
    aybw_term, aybw_term_clamped = np.empty(2), np.empty(2)
    aybw_term[0] = -(3.0*a + 2.0*y + 2.0*b/w)
    aybw_term[1] = -(3.0*a + 2.0*y - 2.0*b/w)
    for i in range(2):
        aybw_term_clamped[i] = aybw_term[i]
        if np.abs(aybw_term_clamped[i]) <= aybw_thresh:
            aybw_term_clamped[i] = 0.0
        elif aybw_term_clamped[i] < 0.0:
            aybw_term_clamped[i] = np.nan

    # Handle no real roots - possibly due to floating point error
    if np.isnan(aybw_term_clamped[0]) and np.isnan(aybw_term_clamped[1]):
        if np.abs(aybw_term[0]) <= np.abs(aybw_term[1]):
            minroot_ind = 0
        else:
            minroot_ind = 1
        thresh_match = np.abs(aybw_term[minroot_ind])/aybw_thresh*threshold_fac
        aybw_term_clamped[minroot_ind] = 0.0
        print(
            "No real roots found in depressed_quartic solver. Smallest term is",
            aybw_term[minroot_ind],
            "(greater than relative threshold for zeroing, which is",
            aybw_thresh,
            "). Bending threshold for this root and forcing to zero.",
            "Threshold is currently", threshold_fac,
            "and would need to be greater than", thresh_match, "to catch this case."
        )

    x1 = -B/(4.0*A) + 0.5*( w + np.sqrt(aybw_term_clamped[0]))
    x2 = -B/(4.0*A) + 0.5*( w - np.sqrt(aybw_term_clamped[0]))
    x3 = -B/(4.0*A) + 0.5*(-w + np.sqrt(aybw_term_clamped[1]))
    x4 = -B/(4.0*A) + 0.5*(-w - np.sqrt(aybw_term_clamped[1]))

    roots = (x1, x2, x3, x4)
    return roots
