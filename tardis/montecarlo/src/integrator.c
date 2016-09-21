#define _USE_MATH_DEFINES

#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "storage.h"
#include "integrator.h"
#include "cmontecarlo.h"

#include <omp.h>

#define NULEN   0
#define LINELEN 1
#define PLEN    2
#define SHELLEN 3

#define C_INV 3.33564e-11
#define M_PI acos (-1)

double integrate_intensity(const double* I_nu, const double h, int N)
{
    double result = (I_nu[0] + I_nu[N-1])/2;
    for (int idx = 1; idx < N-1; ++idx)
    {
        result += I_nu[idx];
    }
    return result*h;
}

void debug_print_arg(double* arg,int len)
{
    for (int64_t i = 0; i < len; i++)
    {
        printf("%e, ",arg[i]);
    }
}
void debug_print_2d_arg(double* arg,int len1, int len2)
{
    for (int64_t i = 0; i < len1; ++i)
    {
        for(int64_t j = 0; j < len2; ++j)
        {
        printf("[%ld,%ld,%ld]: %.8f, ",i,j,i*len2+j,arg[i*len2+j]);
        }
        printf("\n\n");
    }
}

/*

void integrate_source_functions(double* L_nu, const double* line_nu, const double* taus, const double* att_S_ul, const double* I_BB,
        const double* nus, const double* ps, const double* Rs, double R_ph, const int64_t* lens)
{
    double* I_nu  = calloc(lens[PLEN], sizeof(double));
    double  inv_t = 1.8432e+42;
    int no_of_cr_shells;
    indexpair_t nu_lims;
    for (int nu_idx = 0; nu_idx < lens[NULEN]; ++nu_idx)
    {
        memset(I_nu,0.0, lens[PLEN] * sizeof(I_nu));
        for (int p_idx = 0; p_idx < lens[PLEN]; ++p_idx)
        {
            if (ps[p_idx] < R_ph)
            {
                I_nu[p_idx] = I_BB[nu_idx];
            }
            no_of_cr_shells = get_num_shell_cr(ps[p_idx],Rs,lens[SHELLEN]);
            printf("no_sh %d\n", no_of_cr_shells);
            printf("start %d, end %d\n", get_cr_start(no_of_cr_shells, ps[p_idx], R_ph),  2*no_of_cr_shells);
            for (int cr_idx = get_cr_start(no_of_cr_shells, ps[p_idx], R_ph); cr_idx < 2*no_of_cr_shells; ++cr_idx)
            {
                nu_lims = find_nu_limits_for_crossing_and_p(nus[nu_idx], ps[p_idx], cr_idx, no_of_cr_shells, inv_t, Rs, line_nu, lens[LINELEN]);
                printf("Inner loop from %d to %d; ",nu_lims.start,nu_lims.end);
                for (int k_idx = nu_lims.start; k_idx < nu_lims.end; ++k_idx)
                {
                    I_nu[p_idx] = I_nu[p_idx] * exp(-taus[k_idx * lens[SHELLEN] + get_sh_idx(cr_idx,no_of_cr_shells)])
                                    + att_S_ul[k_idx * lens[SHELLEN] + get_sh_idx(cr_idx,no_of_cr_shells)];
                }
            }
        }
        //L_nu[nu_idx] = integrate_intensity(I_nu, ps, lens[PLEN]);
    }
    printf("\n\n");
}
*/

static inline double
calculate_z(double r, double p, double inv_t)
{
    return (r > p) ? sqrt(r * r - p * p) * C_INV * inv_t : 0;
}


int64_t populate_z(const storage_model_t *storage, const double p, double *oz, int64_t *oshell_id)
{
    //const double *radius = storage->r_outer;

    // Abbreviations
    double *r = storage->r_outer;
    const int64_t N = storage->no_of_shells;
    double inv_t = storage->inverse_time_explosion;
    double z = 0;

    int64_t i = 0, offset = -1, i_low, i_up;

    if (p <= storage->r_inner[0])
    {
        oz[0] = calculate_z(storage->r_inner[0], p, inv_t);
        oshell_id[0] = 0;
        for(i = 0; i < N; ++i)
        { // Loop from inside to outside
            oz[i+1] = calculate_z(r[i], p, inv_t);
            oshell_id[i+1] = i;
        }
        return N + 1;
    }
    else
    {
        for(i = 0; i < N; ++i)
        { // Loop from inside to outside
            z = calculate_z(r[i], p, inv_t);
            if (z==0)
                continue;
            if (offset == -1)
            {
                offset = i;
            }
            i_low = N - i - 1;
            i_up = N + i - 2 * offset;

            oz[i_low] = -z;
            oshell_id[i_low] = i;
            oz[i_up] = z;
            oshell_id[i_up] = i;
        }
        return 2*( N - offset);
    }
}


void _formal_integral(
        storage_model_t *storage, double *I_BB, double *att_S_ul, int N, double *L)
{
    // Initialization phase
#pragma omp parallel shared(L)
    {
        printf("Doing the formal integral with %d threads", omp_get_num_threads());

        int64_t offset = 0, i = 0,
                size_line = storage->no_of_lines,
                size_shell = storage->no_of_shells,
                size_tau = size_line * size_shell,
                size_z = 2 * size_shell + 1,
                idx_nu_start = 0;


        double *I_nu  = calloc(N, sizeof(double));
        double *z = calloc(2 * storage->no_of_shells + 1, sizeof(double));
        int64_t *shell_id = calloc(2 * storage->no_of_shells + 1, sizeof(int64_t));
        double *exp_tau = malloc(size_tau * sizeof(double));
        //double exp_tau[size_tau];


        // TODO: This omits the last bin sometimes
        int64_t spectrum_length =
            (storage->spectrum_end_nu - storage->spectrum_start_nu)/storage->spectrum_delta_nu;

        double R_ph = storage->r_inner[0];
        double R_max = storage->r_outer[size_shell - 1];
        double p = 0, nu_start, nu_end, nu, exp_factor;

        double *pexp_tau, *patt_S_ul, *pline;

        // Prepare exp_tau
        for (i = 0; i < size_tau; ++i) {
            exp_tau[i] = exp( -storage->line_lists_tau_sobolevs[i]);
        }

        // Loop over wavelengths in spectrum
#pragma omp for
        for (int nu_idx = 0; nu_idx < spectrum_length ; ++nu_idx)
        {
            nu = storage->spectrum_start_nu + nu_idx * storage->spectrum_delta_nu;

            // Loop over discrete values along line
            for (int p_idx = 0; p_idx < N; ++p_idx)
            {
                // TODO: precompute these and save as 2D array
                memset(z, 0, size_z * sizeof(*z));
                memset(shell_id, 0, size_z * sizeof(*shell_id));

                // Maybe correct? At least this matches the BB *exacly*
                p = R_max/N * (p_idx + 0.5);

                populate_z(storage, p, z, shell_id);

                // initialize I_nu
                if (p <= R_ph)
                    I_nu[p_idx] = I_BB[nu_idx];
                else
                    I_nu[p_idx] = 0;

                // TODO: Ugly loop
                // Loop over all intersections

                // TODO: replace by number of intersections and remove break
                for (i = 0; i < size_z; ++i)
                {
                    if (z[i] == 0)
                        break;
                    nu_start = nu * ( 1 - z[i]);
                    nu_end = nu * ( 1 - z[i+1]);

                    // Calculate offset properly
                    // Which shell is important for photosphere?
                    offset = shell_id[i] * size_line;

                    // Find first contributing line
                    line_search(
                            storage->line_list_nu, nu_start, size_line,
                            &idx_nu_start
                            );

                    // Initialize pointers for inner loop
                    pline = storage->line_list_nu + idx_nu_start;
                    pexp_tau = exp_tau + offset + idx_nu_start;
                    patt_S_ul = att_S_ul + offset + idx_nu_start;

                    for (;pline < storage->line_list_nu + size_line;
                            // We have to increment all pointers simultanously
                            ++pline,
                            ++pexp_tau,
                            ++patt_S_ul)
                    {
                        if (*pline < nu_end)
                            break;
                        I_nu[p_idx] = I_nu[p_idx] * (*pexp_tau) + *patt_S_ul;

                    }
                }
                I_nu[p_idx] *= p;
            }
            L[nu_idx] = 8 * M_PI * M_PI * integrate_intensity(I_nu, R_max/N, N);
        }

        // Free everything allocated on heap
        free(z);
        free(shell_id);
        free(I_nu);
        printf("\n\n");
    }
}
