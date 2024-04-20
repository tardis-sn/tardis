import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.linalg as linalg
from scipy.interpolate import interp1d
from astropy import units as u
from tardis import constants as const
from numba import njit, char, float64, int64, typeof, byte, prange
from numba.experimental import jitclass


from tardis.montecarlo.montecarlo_numba.numba_config import SIGMA_THOMSON
from tardis.montecarlo.montecarlo_numba import njit_dict, njit_dict_no_parallel
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    opacity_state_initialize,
    NumbaModel,
    OpacityState,
)
from tardis.montecarlo.montecarlo_numba.formal_integral_cuda import (
    CudaFormalIntegrator,
)

from tardis.montecarlo.spectrum import TARDISSpectrum

C_INV = 3.33564e-11
M_PI = np.arccos(-1)
KB_CGS = 1.3806488e-16
H_CGS = 6.62606957e-27


class IntegrationError(Exception):
    pass


@njit(**njit_dict)
def numba_formal_integral(
    geometry,
    model,
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
    model, plasma, and estimator are the numba variants

    Returns
    -------
    L : float64 array
        integrated luminosities
    I_nu_p : float64 2D array
        intensities at each p-ray multiplied by p
        frequency x p-ray grid
    """

    # todo: add all the original todos
    # Initialize the output which is shared among threads
    L = np.zeros(inu_size, dtype=np.float64)
    # global read-only values
    size_line, size_shell = tau_sobolev.shape
    size_tau = size_line * size_shell
    R_ph = geometry.r_inner[0]  # make sure these are cgs
    R_max = geometry.r_outer[size_shell - 1]
    pp = np.zeros(N, dtype=np.float64)  # check
    exp_tau = np.zeros(size_tau, dtype=np.float64)
    exp_tau = np.exp(-tau_sobolev.T.ravel())  # maybe make this 2D?
    pp[::] = calculate_p_values(R_max, N)
    line_list_nu = plasma.line_list_nu
    # done with instantiation
    # now loop over wavelength in spectrum
    I_nu_p = np.zeros((inu_size, N), dtype=np.float64)
    for nu_idx in prange(inu_size):
        I_nu = I_nu_p[nu_idx]
        z = np.zeros(2 * size_shell, dtype=np.float64)
        shell_id = np.zeros(2 * size_shell, dtype=np.int64)
        offset = 0
        size_z = 0
        idx_nu_start = 0
        direction = 0
        first = 0
        i = 0
        p = 0.0
        nu_start = 0.0
        nu_end = 0.0
        nu = 0.0
        zstart = 0.0
        zend = 0.0
        escat_contrib = 0.0
        escat_op = 0.0
        Jkkp = 0.0
        pexp_tau = 0
        patt_S_ul = 0
        pJred_lu = 0
        pJblue_lu = 0
        pline = 0

        nu = inu[nu_idx]
        # now loop over discrete values along line
        for p_idx in range(1, N):
            escat_contrib = 0
            p = pp[p_idx]

            # initialize z intersections for p values
            size_z = populate_z(
                geometry, model, p, z, shell_id
            )  # check returns
            # initialize I_nu
            if p <= R_ph:
                I_nu[p_idx] = intensity_black_body(nu * z[0], iT)
            else:
                I_nu[p_idx] = 0

            # find first contributing lines
            nu_start = nu * z[0]
            nu_end = nu * z[1]
            idx_nu_start = line_search(plasma.line_list_nu, nu_start, size_line)
            offset = shell_id[0] * size_line
            # start tracking accumulated e-scattering optical depth
            zstart = model.time_explosion / C_INV * (1.0 - z[0])
            # Initialize "pointers"
            pline = int(idx_nu_start)
            pexp_tau = int(offset + idx_nu_start)
            patt_S_ul = int(offset + idx_nu_start)
            pJred_lu = int(offset + idx_nu_start)
            pJblue_lu = int(offset + idx_nu_start)

            # flag for first contribution to integration on current p-ray
            first = 1
            nu_ends = nu * z[1:]
            nu_ends_idxs = size_line - np.searchsorted(
                line_list_nu[::-1], nu_ends, side="right"
            )
            # loop over all interactions
            for i in range(size_z - 1):
                escat_op = electron_density[int(shell_id[i])] * SIGMA_THOMSON
                nu_end = nu_ends[i]
                nu_end_idx = nu_ends_idxs[i]
                for _ in range(max(nu_end_idx - pline, 0)):
                    # calculate e-scattering optical depth to next resonance point
                    zend = (
                        model.time_explosion
                        / C_INV
                        * (1.0 - line_list_nu[pline] / nu)
                    )  # check

                    if first == 1:
                        # first contribution to integration
                        # NOTE: this treatment of I_nu_b (given
                        #   by boundary conditions) is not in Lucy 1999;
                        #   should be re-examined carefully
                        escat_contrib += (
                            (zend - zstart)
                            * escat_op
                            * (Jblue_lu[pJblue_lu] - I_nu[p_idx])
                        )
                        first = 0
                    else:
                        # Account for e-scattering, c.f. Eqs 27, 28 in Lucy 1999
                        Jkkp = 0.5 * (Jred_lu[pJred_lu] + Jblue_lu[pJblue_lu])
                        escat_contrib += (
                            (zend - zstart) * escat_op * (Jkkp - I_nu[p_idx])
                        )
                        # this introduces the necessary ffset of one element between
                        # pJblue_lu and pJred_lu
                        pJred_lu += 1
                    I_nu[p_idx] += escat_contrib
                    # // Lucy 1999, Eq 26
                    I_nu[p_idx] *= exp_tau[pexp_tau]
                    I_nu[p_idx] += att_S_ul[patt_S_ul]

                    # // reset e-scattering opacity
                    escat_contrib = 0
                    zstart = zend

                    pline += 1
                    pexp_tau += 1
                    patt_S_ul += 1
                    pJblue_lu += 1

                # calculate e-scattering optical depth to grid cell boundary

                Jkkp = 0.5 * (Jred_lu[pJred_lu] + Jblue_lu[pJblue_lu])
                zend = (
                    model.time_explosion / C_INV * (1.0 - nu_end / nu)
                )  # check
                escat_contrib += (
                    (zend - zstart) * escat_op * (Jkkp - I_nu[p_idx])
                )
                zstart = zend

                # advance pointers
                direction = int((shell_id[i + 1] - shell_id[i]) * size_line)
                pexp_tau += direction
                patt_S_ul += direction
                pJred_lu += direction
                pJblue_lu += direction
            I_nu[p_idx] *= p
        L[nu_idx] = 8 * M_PI * M_PI * trapezoid_integration(I_nu, R_max / N)

    return L, I_nu_p


# integrator_spec = [
#    ("model", NumbaModel.class_type.instance_type),
#    ("plasma", OpacityState.class_type.instance_type),
#    ("points", int64),
# ]


# @jitclass(integrator_spec)
class NumbaFormalIntegrator(object):
    """
    Helper class for performing the formal integral
    with numba.
    """

    def __init__(self, geometry, model, plasma, points=1000):
        self.geometry = geometry
        self.model = model
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
            self.model,
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


class FormalIntegrator(object):
    """
    Class containing the formal integrator.

    If there is a NVIDIA CUDA GPU available,
    the formal integral will automatically run
    on it. If multiple GPUs are available, it will
    choose the first one that it sees. You can
    read more about selecting different GPUs on
    Numba's CUDA documentation.

    Parameters
    ----------
    model : tardis.model.SimulationState
    plasma : tardis.plasma.BasePlasma
    transport : tardis.montecarlo.MontecarloTransport
    points : int64
    """

    def __init__(self, simulation_state, plasma, transport, points=1000):
        self.simulation_state = simulation_state
        self.transport = transport
        self.points = points
        if transport:
            self.montecarlo_configuration = (
                self.transport.montecarlo_configuration
            )
        if plasma:
            self.plasma = opacity_state_initialize(
                plasma,
                transport.line_interaction_type,
                self.montecarlo_configuration.DISABLE_LINE_SCATTERING,
                self.montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED,
            )
            self.atomic_data = plasma.atomic_data
            self.original_plasma = plasma
            self.levels_index = plasma.levels

    def generate_numba_objects(self):
        """instantiate the numba interface objects
        needed for computing the formal integral"""
        from tardis.model.geometry.radial1d import NumbaRadial1DGeometry

        self.numba_radial_1d_geometry = NumbaRadial1DGeometry(
            self.transport.r_inner_i,
            self.transport.r_outer_i,
            self.transport.r_inner_i
            / self.simulation_state.time_explosion.to("s").value,
            self.transport.r_outer_i
            / self.simulation_state.time_explosion.to("s").value,
        )
        self.numba_model = NumbaModel(
            self.simulation_state.time_explosion.cgs.value,
        )
        self.opacity_state = opacity_state_initialize(
            self.original_plasma,
            self.transport.line_interaction_type,
            self.montecarlo_configuration.DISABLE_LINE_SCATTERING,
            self.montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED,
        )
        if self.transport.use_gpu:
            self.integrator = CudaFormalIntegrator(
                self.numba_radial_1d_geometry,
                self.numba_model,
                self.opacity_state,
                self.points,
            )
        else:
            self.integrator = NumbaFormalIntegrator(
                self.numba_radial_1d_geometry,
                self.numba_model,
                self.opacity_state,
                self.points,
            )

    def check(self, raises=True):
        """
        A method that determines if the formal integral can be performed with
        the current configuration settings

        The function returns False if the configuration conflicts with the
        required settings. If raises evaluates to True, then a
        IntegrationError is raised instead
        """

        def raise_or_return(message):
            if raises:
                raise IntegrationError(message)
            else:
                warnings.warn(message)
                return False

        for obj in (self.simulation_state, self.plasma, self.transport):
            if obj is None:
                return raise_or_return(
                    "The integrator is missing either model, plasma or "
                    "transport. Please make sure these are provided to the "
                    "FormalIntegrator."
                )

        if not self.transport.line_interaction_type in [
            "downbranch",
            "macroatom",
        ]:
            return raise_or_return(
                "The FormalIntegrator currently only works for "
                'line_interaction_type == "downbranch"'
                'and line_interaction_type == "macroatom"'
            )

        if self.montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED:
            return raise_or_return(
                "The FormalIntegrator currently does not work for "
                "continuum interactions."
            )

        return True

    def calculate_spectrum(
        self, frequency, points=None, interpolate_shells=0, raises=True
    ):
        # Very crude implementation
        # The c extension needs bin centers (or something similar)
        # while TARDISSpectrum needs bin edges
        self.check(raises)
        N = points or self.points
        if interpolate_shells == 0:  # Default Value
            interpolate_shells = max(2 * self.simulation_state.no_of_shells, 80)
            warnings.warn(
                "The number of interpolate_shells was not "
                f"specified. The value was set to {interpolate_shells}."
            )
        self.interpolate_shells = interpolate_shells
        frequency = frequency.to("Hz", u.spectral())

        luminosity = u.Quantity(self.formal_integral(frequency, N), "erg") * (
            frequency[1] - frequency[0]
        )

        # Ugly hack to convert to 'bin edges'
        frequency = u.Quantity(
            np.concatenate(
                [
                    frequency.value,
                    [frequency.value[-1] + np.diff(frequency.value)[-1]],
                ]
            ),
            frequency.unit,
        )

        return TARDISSpectrum(frequency, luminosity)

    def make_source_function(self):
        """
        Calculates the source function using the line absorption rate estimator `Edotlu_estimator`

        Formally it calculates the expression ( 1 - exp(-tau_ul) ) S_ul but this product is what we need later,
        so there is no need to factor out the source function explicitly.

        Parameters
        ----------
        model : tardis.model.SimulationState

        Returns
        -------
        Numpy array containing ( 1 - exp(-tau_ul) ) S_ul ordered by wavelength of the transition u -> l
        """
        simulation_state = self.simulation_state
        transport = self.transport
        montecarlo_transport_state = transport.transport_state

        # macro_ref = self.atomic_data.macro_atom_references
        macro_ref = self.atomic_data.macro_atom_references
        # macro_data = self.atomic_data.macro_atom_data
        macro_data = self.original_plasma.macro_atom_data

        no_lvls = len(self.levels_index)
        no_shells = len(simulation_state.dilution_factor)

        if transport.line_interaction_type == "macroatom":
            internal_jump_mask = (macro_data.transition_type >= 0).values
            ma_int_data = macro_data[internal_jump_mask]
            internal = self.original_plasma.transition_probabilities[
                internal_jump_mask
            ]

            source_level_idx = ma_int_data.source_level_idx.values
            destination_level_idx = ma_int_data.destination_level_idx.values

        Edotlu_norm_factor = 1 / (
            montecarlo_transport_state.packet_collection.time_of_simulation
            * simulation_state.volume
        )
        exptau = 1 - np.exp(-self.original_plasma.tau_sobolevs)
        Edotlu = (
            Edotlu_norm_factor
            * exptau
            * montecarlo_transport_state.radfield_mc_estimators.Edotlu_estimator
        )

        # The following may be achieved by calling the appropriate plasma
        # functions
        Jbluelu_norm_factor = (
            (
                const.c.cgs
                * simulation_state.time_explosion
                / (
                    4
                    * np.pi
                    * montecarlo_transport_state.time_of_simulation
                    * simulation_state.volume
                )
            )
            .to("1/(cm^2 s)")
            .value
        )
        # Jbluelu should already by in the correct order, i.e. by wavelength of
        # the transition l->u
        Jbluelu = (
            transport.transport_state.radfield_mc_estimators.j_blue_estimator
            * Jbluelu_norm_factor
        )

        upper_level_index = self.atomic_data.lines.index.droplevel(
            "level_number_lower"
        )
        e_dot_lu = pd.DataFrame(Edotlu.values, index=upper_level_index)
        e_dot_u = e_dot_lu.groupby(level=[0, 1, 2]).sum()
        e_dot_u_src_idx = macro_ref.loc[e_dot_u.index].references_idx.values

        if transport.line_interaction_type == "macroatom":
            C_frame = pd.DataFrame(
                columns=np.arange(no_shells), index=macro_ref.index
            )
            q_indices = (source_level_idx, destination_level_idx)
            for shell in range(no_shells):
                Q = sp.coo_matrix(
                    (internal[shell], q_indices), shape=(no_lvls, no_lvls)
                )
                inv_N = sp.identity(no_lvls) - Q
                e_dot_u_vec = np.zeros(no_lvls)
                e_dot_u_vec[e_dot_u_src_idx] = e_dot_u[shell].values
                C_frame[shell] = linalg.spsolve(inv_N.T, e_dot_u_vec)

        e_dot_u.index.names = [
            "atomic_number",
            "ion_number",
            "source_level_number",
        ]  # To make the q_ul e_dot_u product work, could be cleaner
        transitions = self.original_plasma.atomic_data.macro_atom_data[
            self.original_plasma.atomic_data.macro_atom_data.transition_type
            == -1
        ].copy()
        transitions_index = transitions.set_index(
            ["atomic_number", "ion_number", "source_level_number"]
        ).index.copy()
        tmp = self.original_plasma.transition_probabilities[
            (self.atomic_data.macro_atom_data.transition_type == -1).values
        ]
        q_ul = tmp.set_index(transitions_index)
        t = simulation_state.time_explosion.value
        t = simulation_state.time_explosion.value
        lines = self.atomic_data.lines.set_index("line_id")
        wave = lines.wavelength_cm.loc[
            transitions.transition_line_id
        ].values.reshape(-1, 1)
        if transport.line_interaction_type == "macroatom":
            e_dot_u = C_frame.loc[e_dot_u.index]
        att_S_ul = wave * (q_ul * e_dot_u) * t / (4 * np.pi)

        result = pd.DataFrame(
            att_S_ul.values, index=transitions.transition_line_id.values
        )
        att_S_ul = result.loc[lines.index.values].values

        # Jredlu should already by in the correct order, i.e. by wavelength of
        # the transition l->u (similar to Jbluelu)
        Jredlu = Jbluelu * np.exp(-self.original_plasma.tau_sobolevs) + att_S_ul
        if self.interpolate_shells > 0:
            (
                att_S_ul,
                Jredlu,
                Jbluelu,
                e_dot_u,
            ) = self.interpolate_integrator_quantities(
                att_S_ul, Jredlu, Jbluelu, e_dot_u
            )
        else:
            transport.r_inner_i = (
                montecarlo_transport_state.geometry_state.r_inner
            )
            transport.r_outer_i = (
                montecarlo_transport_state.geometry_state.r_outer
            )
            transport.tau_sobolevs_integ = (
                self.original_plasma.tau_sobolevs.values
            )
            transport.electron_densities_integ = (
                self.original_plasma.electron_densities.values
            )

        return att_S_ul, Jredlu, Jbluelu, e_dot_u

    def interpolate_integrator_quantities(
        self, att_S_ul, Jredlu, Jbluelu, e_dot_u
    ):
        transport = self.transport
        mct_state = transport.transport_state
        plasma = self.original_plasma
        nshells = self.interpolate_shells
        r_middle = (
            mct_state.geometry_state.r_inner + mct_state.geometry_state.r_outer
        ) / 2.0

        r_integ = np.linspace(
            mct_state.geometry_state.r_inner[0],
            mct_state.geometry_state.r_outer[-1],
            nshells,
        )
        transport.r_inner_i = r_integ[:-1]
        transport.r_outer_i = r_integ[1:]

        r_middle_integ = (r_integ[:-1] + r_integ[1:]) / 2.0

        transport.electron_densities_integ = interp1d(
            r_middle,
            plasma.electron_densities,
            fill_value="extrapolate",
            kind="nearest",
        )(r_middle_integ)
        # Assume tau_sobolevs to be constant within a shell
        # (as in the MC simulation)
        transport.tau_sobolevs_integ = interp1d(
            r_middle,
            plasma.tau_sobolevs,
            fill_value="extrapolate",
            kind="nearest",
        )(r_middle_integ)
        att_S_ul = interp1d(r_middle, att_S_ul, fill_value="extrapolate")(
            r_middle_integ
        )
        Jredlu = pd.DataFrame(
            interp1d(r_middle, Jredlu, fill_value="extrapolate")(r_middle_integ)
        )
        Jbluelu = interp1d(r_middle, Jbluelu, fill_value="extrapolate")(
            r_middle_integ
        )
        e_dot_u = interp1d(r_middle, e_dot_u, fill_value="extrapolate")(
            r_middle_integ
        )

        # Set negative values from the extrapolation to zero
        att_S_ul = att_S_ul.clip(0.0)
        Jbluelu = Jbluelu.clip(0.0)
        Jredlu = Jredlu.clip(0.0)
        e_dot_u = e_dot_u.clip(0.0)
        return att_S_ul, Jredlu, Jbluelu, e_dot_u

    def formal_integral(self, nu, N):
        """Do the formal integral with the numba
        routines"""
        # TODO: get rid of storage later on

        res = self.make_source_function()

        att_S_ul = res[0].flatten(order="F")
        Jred_lu = res[1].values.flatten(order="F")
        Jblue_lu = res[2].flatten(order="F")

        self.generate_numba_objects()
        L, I_nu_p = self.integrator.formal_integral(
            self.simulation_state.t_inner,
            nu,
            nu.shape[0],
            att_S_ul,
            Jred_lu,
            Jblue_lu,
            self.transport.tau_sobolevs_integ,
            self.transport.electron_densities_integ,
            N,
        )
        R_max = self.transport.r_outer_i[-1]
        ps = calculate_p_values(R_max, N)[None, :]
        I_nu_p[:, 1:] /= ps[:, 1:]
        self.transport.I_nu_p = I_nu_p
        self.transport.p_rays = ps

        I_nu = self.transport.I_nu_p * ps
        L_test = np.array(
            [
                8 * M_PI * M_PI * trapezoid_integration((I_nu)[i, :], R_max / N)
                for i in range(nu.shape[0])
            ]
        )
        error = np.max(np.abs((L_test - L) / L))
        assert (
            error < 1e-7
        ), f"Incorrect I_nu_p values, max relative difference:{error}"

        return np.array(L, np.float64)


@njit(**njit_dict_no_parallel)
def populate_z(geometry, model, p, oz, oshell_id):
    """Calculate p line intersections

    This function calculates the intersection points of the p-line with
    each shell

    Inputs:
        :p: (double) distance of the integration line to the center
        :oz: (array of doubles) will be set with z values. the array is truncated
                    by the value `1`.
        :oshell_id: (int64) will be set with the corresponding shell_ids
    """
    # abbreviations
    r = geometry.r_outer
    N = len(geometry.r_inner)  # check
    inv_t = 1 / model.time_explosion
    z = 0
    offset = N

    if p <= geometry.r_inner[0]:
        # intersect the photosphere
        for i in range(N):
            oz[i] = 1 - calculate_z(r[i], p, inv_t)
            oshell_id[i] = i
        return N
    else:
        # no intersection with photosphere
        # that means we intersect each shell twice
        for i in range(N):
            z = calculate_z(r[i], p, inv_t)
            if z == 0:
                continue
            if offset == N:
                offset = i
            # calculate the index in the resulting array
            i_low = N - i - 1  # the far intersection with the shell
            i_up = N + i - 2 * offset  # the nearer intersection with the shell

            # setting the arrays; check return them?
            oz[i_low] = 1 + z
            oshell_id[i_low] = i
            oz[i_up] = 1 - z
            oshell_id[i_up] = i
        return 2 * (N - offset)


@njit(**njit_dict_no_parallel)
def calculate_z(r, p, inv_t):
    """Calculate distance to p line

    Calculate half of the length of the p-line inside a shell
    of radius r in terms of unit length (c * t_exp).
    If shell and p-line do not intersect, return 0.

    Inputs:
        :r: (double) radius of the shell
        :p: (double) distance of the p-line to the center of the supernova
        :inv_t: (double) inverse time_explosio is needed to norm to unit-length
    """
    if r > p:
        return np.sqrt(r * r - p * p) * C_INV * inv_t
    else:
        return 0


class BoundsError(IndexError):
    pass


@njit(**njit_dict_no_parallel)
def line_search(nu, nu_insert, number_of_lines):
    """
    Insert a value in to an array of line frequencies

    Inputs:
        :nu: (array) line frequencies
        :nu_insert: (int) value of nu key
        :number_of_lines: (int) number of lines in the line list

    Outputs:
        index of the next line ot the red.
                If the key value is redder
                 than the reddest line returns number_of_lines.
    """
    # TODO: fix the TARDIS_ERROR_OK
    # tardis_error_t ret_val = TARDIS_ERROR_OK # check
    imin = 0
    imax = number_of_lines - 1
    if nu_insert > nu[imin]:
        result = imin
    elif nu_insert < nu[imax]:
        result = imax + 1
    else:
        result = reverse_binary_search(nu, nu_insert, imin, imax)
        result = result + 1
    return result


@njit(**njit_dict_no_parallel)
def reverse_binary_search(x, x_insert, imin, imax):
    """Look for a place to insert a value in an inversely sorted float array.

    Inputs:
        :x: (array) an inversely (largest to lowest) sorted float array
        :x_insert: (value) a value to insert
        :imin: (int) lower bound
        :imax: (int) upper bound

    Outputs:
        index of the next boundary to the left
    """
    # ret_val = TARDIS_ERROR_OK # check
    if (x_insert > x[imin]) or (x_insert < x[imax]):
        raise BoundsError  # check
    return len(x) - 1 - np.searchsorted(x[::-1], x_insert, side="right")


@njit(**njit_dict_no_parallel)
def trapezoid_integration(array, h):
    """in the future, let's just replace
    this with the numpy trapz
    since it is numba compatable
    """
    return np.trapz(array, dx=h)


@njit(**njit_dict_no_parallel)
def intensity_black_body(nu, temperature):
    """
    Calculate the blackbody intensity.

    Parameters
    ----------
    nu : float64
        frequency
    temperature : float64
        Temperature

    Returns
    -------
    float64
    """
    if nu == 0:
        return np.nan  # to avoid ZeroDivisionError
    beta_rad = 1 / (KB_CGS * temperature)
    coefficient = 2 * H_CGS * C_INV * C_INV
    return coefficient * nu * nu * nu / (np.exp(H_CGS * nu * beta_rad) - 1)


@njit(**njit_dict_no_parallel)
def calculate_p_values(R_max, N):
    """This can probably be replaced with a simpler function"""
    return np.arange(N).astype(np.float64) * R_max / (N - 1)
