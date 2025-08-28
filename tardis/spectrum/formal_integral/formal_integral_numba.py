import numpy as np
from numba import njit, prange

from tardis.spectrum.formal_integral.base import (
    C_INV,
    BoundsError,
    trapezoid_integration,
    setup_formal_integral_inputs
)
from tardis.transport.montecarlo import njit_dict, njit_dict_no_parallel
from tardis.transport.montecarlo.configuration.constants import SIGMA_THOMSON

@njit(**njit_dict_no_parallel)
def increment_escat_contrib(
    escat_contrib,
    first,
    pJred_lu,
    zend,
    zstart,
    escat_op,
    Jblue_lu,
    Jred_lu,
    I_nu_p,
):

    if first == 1:
        # first contribution to integration
        # NOTE: this treatment of I_nu_b (given
        #   by boundary conditions) is not in Lucy 1999;
        #   should be re-examined carefully
        escat_contrib += (zend - zstart) * escat_op * (Jblue_lu - I_nu_p)
        first = 0
    else:
        # Account for e-scattering, c.f. Eqs 27, 28 in Lucy 1999
        Jkkp = 0.5 * (Jred_lu + Jblue_lu)
        escat_contrib += (zend - zstart) * escat_op * (Jkkp - I_nu_p)
        # this introduces the necessary ffset of one element between
        # pJblue_lu and pJred_lu
        pJred_lu += 1

    return escat_contrib, first, pJred_lu

@njit(**njit_dict)
def numba_formal_integral(
    geometry,
    time_explosion,
    plasma,
    iT,
    inu,
    inu_size,
    att_S_ul,
    Jred_lu,
    Jblue_lu,
    tau_sobolev,
    electron_density,
    N,
):
    """
    Returns
    -------
    L : float64 array
        integrated luminosities
    I_nu_p : float64 2D array
        intensities at each p-ray multiplied by p
        frequency x p-ray grid
    """
    # Initialize the output which is shared among threads
    L = np.zeros(inu_size, dtype=np.float64)
    # global read-only values
    # size_line, size_shell = tau_sobolev.shape
    # size_tau = size_line * size_shell
    # R_ph = geometry.r_inner[0]  # make sure these are cgs
    # R_max = geometry.r_outer[size_shell - 1]
    # pp = np.zeros(N, dtype=np.float64)  # check
    # exp_tau = np.zeros(size_tau, dtype=np.float64)
    # exp_tau = np.exp(-tau_sobolev.T.ravel())  # maybe make this 2D?
    # pp[::] = calculate_p_values(R_max, N)
    line_list_nu = plasma.line_list_nu
    
    # I_nu_p, zs, shell_ids, size_zs = init_Inup_numba(
    #     inu, pp, iT, R_ph, N, size_shell, geometry, time_explosion
    # )
    R_max = geometry.r_outer[-1]
    (I_nu_p, ps, size_zs, lines_idx, lines_idx_offset, nu_ends, nu_ends_idxs, zstarts, escat_ops, directions, exp_tau ) = setup_formal_integral_inputs(
            inu, iT, N, geometry, time_explosion, line_list_nu, tau_sobolev, electron_density)


    # done with instantiation
    # now loop over wavelength in spectrum
    for nu_idx in prange(inu_size):
        I_nu = I_nu_p[nu_idx]
        nu = inu[nu_idx]

        # now loop over discrete values along line
        for p_idx in range(1, N):
            escat_contrib = 0
            p = ps[p_idx]
            # z = zs[p_idx]
            # shell_id = shell_ids[p_idx]
            size_z = size_zs[p_idx]

            # find first contributing lines
            # nu_start = nu * z[0]
            # idx_nu_start = line_search(plasma.line_list_nu, nu_start, size_line)
            # offset = shell_id[0] * size_line
            # start tracking accumulated e-scattering optical depth
            # zstart = time_explosion / C_INV * (1.0 - z[0])
            zstart = zstarts[p_idx]

            # Initialize "pointers"
            # pline = int(idx_nu_start)
            # pline_offset = int(offset + idx_nu_start)
            # pJred_lu = int(offset + idx_nu_start)
            line_idx = lines_idx[nu_idx, p_idx]
            line_idx_offset = lines_idx_offset[nu_idx, p_idx]
            line_Jred_lu_idx = lines_idx_offset[nu_idx, p_idx]


            # flag for first contribution to integration on current p-ray
            first = 1
            # nu_ends = nu * z[1:]
            # nu_ends = nu_ends
            # nu_ends_idxs = size_line - np.searchsorted(
            #     line_list_nu[::-1], nu_ends, side="right"
            # )
            # nu_ends_idxs = nu_ends_idxs[nu_idx, p_idx]

            # loop over all interactions
            for i in range(size_z - 1):
                # escat_op = electron_density[int(shell_id[i])] * SIGMA_THOMSON
                escat_op = escat_ops[nu_idx, p_idx, i]
                nu_end = nu_ends[nu_idx, p_idx, i]
                nu_end_idx = nu_ends_idxs[nu_idx, p_idx, i]
                for _ in range(max(nu_end_idx - line_idx, 0)):
                    # calculate e-scattering optical depth to next resonance point
                    zend = (
                        time_explosion
                        / C_INV
                        * (1.0 - line_list_nu[line_idx] / nu)
                    )  # check

                    try:
                        Jblue_lu[line_idx_offset]
                    except:
                        print("line_idx_offset", line_idx_offset)

                    escat_contrib, first, line_Jred_lu_idx = increment_escat_contrib(
                        escat_contrib,
                        first,
                        line_Jred_lu_idx,
                        zend,
                        zstart,
                        escat_op,
                        Jblue_lu[line_idx_offset],
                        Jred_lu[line_Jred_lu_idx],
                        I_nu[p_idx],
                    )

                    I_nu[p_idx] += escat_contrib
                    # // Lucy 1999, Eq 26
                    I_nu[p_idx] *= exp_tau[line_idx_offset]
                    I_nu[p_idx] += att_S_ul[line_idx_offset]

                    # // reset e-scattering opacity
                    escat_contrib = 0
                    zstart = zend

                    line_idx += 1
                    line_idx_offset += 1

                # calculate e-scattering optical depth to grid cell boundary
                Jkkp = 0.5 * (Jred_lu[line_Jred_lu_idx] + Jblue_lu[line_idx_offset])
                zend = time_explosion / C_INV * (1.0 - nu_end / nu)  # check
                escat_contrib += (
                    (zend - zstart) * escat_op * (Jkkp - I_nu[p_idx])
                ) # check if cant use increment function here
                zstart = zend

                # advance pointers
                direction = directions[nu_idx, p_idx, i]
                line_idx_offset += direction
                line_Jred_lu_idx += direction
            I_nu[p_idx] *= p
        L[nu_idx] = 8 * np.pi * np.pi * trapezoid_integration(I_nu, R_max / N)

    return L, I_nu_p


class NumbaFormalIntegrator:
    """
    Helper class for performing the formal integral
    with numba.
    """

    def __init__(self, geometry, time_explosion, plasma, points=1000):
        self.geometry = geometry
        self.time_explosion = time_explosion
        self.plasma = plasma
        self.points = points

    def formal_integral(
        self,
        iT,
        inu,
        inu_size,
        att_S_ul,
        Jred_lu,
        Jblue_lu,
        tau_sobolev,
        electron_density,
        N,
    ):
        """
        Simple wrapper for the numba implementation of the formal integral
        """
        return numba_formal_integral(
            self.geometry,
            self.time_explosion,
            self.plasma,
            iT,
            inu,
            inu_size,
            att_S_ul,
            Jred_lu,
            Jblue_lu,
            tau_sobolev,
            electron_density,
            N,
        )
