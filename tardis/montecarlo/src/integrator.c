#include <string.h>
#include <stdio.h>
#include <math.h>

#include "integrator.h"

#ifdef WITHOPENMP
#include <omp.h>
#endif

#define NULEN   0
#define LINELEN 1
#define PLEN    2
#define SHELLEN 3

indexpair_t find_array_bounds_nonzero(double* array, int idx2, int len)
{
    indexpair_t idxs;
    idxs.start = -1;
    idxs.end   = -1;
    for (int i=0; i < len; ++i){
        if (idxs.start != -1){
            if (array[i,idx2] > 0){
                idxs.start = i; 
                idxs.end   = i;
            }
        }
        if (array[i,idx2] > 0){
             idxs.end   = i;
        }
    }
    return idxs;
}

indexpair_t find_nu_limits_for_shell_and_p(double p,int shl_idx, double* Rs, int no_shells)
{
    assert(Rs[0] > Rs[1]); // Decreasing order
    double R_ph = Rs[no_shells];
    printf("Not implemented yet");
}

double integrate_intensity(double* I_nu,double* ps, int* lens)
{
    printf("Not implemented yet");
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
        printf("[%d,%d,%d]: %.8f, ",i,j,i*len2+j,arg[i*len2+j]);
        }
        printf("\n\n");
    }
}


void integrate_source_functions(double* L_nu, double* line_nu, double* taus, double* att_S_ul, double* I_BB, double* nus, 
              double* ps, double* Rs, double R_ph, double ct, int* lens)
{
    double I_nu[lens[SHELLEN]];
    indexpair_t nu_lims;
    for (int nu_idx = 0; nu_idx < lens[NULEN]; ++nu_idx)
    {
        memset(I_nu,0.0,sizeof(I_nu));
        for (int p_idx = 0; p_idx < lens[PLEN]; ++p_idx)
        {
            if (ps[p_idx] < R_ph) 
            {
                I_nu[nu_idx] = I_BB[nu_idx];
            }
            for (int shl_idx = 0; shl_idx < lens[SHELLEN]; ++shl_idx)
            {
                nu_lims = find_nu_limits_for_shell_and_p(ps[p_idx],shl_idx,Rs,R_ph);
                for (int k_idx = nu_lims.start; k_idx < nu_lims.end; ++k_idx)
                {
                    I_nu[nu_idx] = I_nu[nu_idx]*exp(-taus[k_idx*lens[SHELLEN]+shl_idx]) + att_S_ul[k_idx*lens[SHELLEN]+shl_idx];
                }
            }
        }
        L_nu[nu_idx] = integrate_intensity(I_nu, ps, lens); 
    }
}
