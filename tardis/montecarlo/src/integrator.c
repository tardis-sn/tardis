#include <string.h>
#include <assert.h>
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

#define M_PI acos(-1.0)
#define C_INV 3.33564e-11

indexpair_t find_nu_limits_for_crossing_and_p(double nu, double p, int cr_idx, int no_of_cr_shells, double inv_ct, const double* Rs, const double* line_nu, int len)
{
    double blu_R, red_R, z_blu, z_red, z_cr, nu_blu, nu_red;
    indexpair_t pair;

    if (no_of_cr_shells > 1)
    {
        assert(Rs[0] > Rs[1]); // Decreasing order
        
        blu_R = get_r(cr_idx,no_of_cr_shells,Rs);
        red_R = get_r(cr_idx+1,no_of_cr_shells,Rs);
        z_blu = sqrt( blu_R*blu_R - p*p );
        z_red = sqrt( red_R*red_R - p*p );
        nu_blu = nu * (1 + get_cr_sign(cr_idx,no_of_cr_shells)*z_blu*inv_ct);
        nu_red = nu * (1 + get_cr_sign(cr_idx,no_of_cr_shells)*z_red*inv_ct);
    }
    else 
    {
        z_cr = sqrt( Rs[cr_idx]*Rs[cr_idx] - p*p );
        nu_blu = nu * (1 - z_cr*inv_ct);
        nu_red = nu * (1 + z_cr*inv_ct);
    }

    for (int idx = 0; idx < len; ++idx)
    {
        if (line_nu[idx] >= nu_blu){
            pair.start = idx;}
        if (line_nu[idx] <= nu_red){
            pair.end = idx-1;}
    }
    
    return pair;
}

double test_nu_limits_for_crossing_and_p(double nu, double p, int cr_idx, int no_of_cr_shells, double inv_ct, const double* Rs, const double* line_nu, int len)
{
    double blu_R, red_R, z_blu, z_red, z_cr, nu_blu, nu_red;
    double pair;

    if (no_of_cr_shells > 1)
    {
        assert(Rs[0] > Rs[1]); // Decreasing order
        
        blu_R = get_r(cr_idx,no_of_cr_shells,Rs);
        red_R = get_r(cr_idx+1,no_of_cr_shells,Rs);
        z_blu = sqrt( blu_R*blu_R - p*p );
        z_red = sqrt( red_R*red_R - p*p );
        nu_blu = nu * (1 + get_cr_sign(cr_idx,no_of_cr_shells)*z_blu*inv_ct);
        nu_red = nu * (1 + get_cr_sign(cr_idx,no_of_cr_shells)*z_red*inv_ct);
    }
    else 
    {
        z_cr = sqrt( Rs[cr_idx]*Rs[cr_idx] - p*p );
        nu_blu = nu * (1 - z_cr*inv_ct);
        nu_red = nu * (1 + z_cr*inv_ct);
    }
    pair =  (1 - z_blu*inv_ct);
    
    return nu_blu;
}


double get_r(int cr_idx, int no_of_cr_shells, const double* Rs)
{
    return Rs[ get_sh_idx(cr_idx, no_of_cr_shells)];
}

int get_cr_sign(int cr_idx, const int no_of_cr_shells)
{
    if (cr_idx < no_of_cr_shells){
        return -1;}
    else if (cr_idx < 2*no_of_cr_shells){
        return 1;}
}

int get_cr_start(int no_of_cr_shells, double p, double R_ph)
{
    if (p >= R_ph) {
        return 0;}
    else if (p < R_ph) {
        return no_of_cr_shells;
    }
}

int get_sh_idx(int cr_idx, int no_of_cr_shells)
{
    if (cr_idx < no_of_cr_shells){
       return cr_idx;}
    else if (cr_idx == 2*no_of_cr_shells-1){ // last crossing, used for indexing Rs
        return 0;}
    else if (cr_idx < 2*no_of_cr_shells-1){
        return (2*(no_of_cr_shells-1) - cr_idx);}

}
int get_num_shell_cr(double p, const double* Rs, int len)
{
    assert(Rs[0] > Rs[1]);

    int num;
    for(num = 0; num < len; ++num)
    {
        if( p >= Rs[num]) {
            break;}
    }
    return num;
}


double integrate_intensity(const double* I_nu, const double* ps, int len)
{
    double result = 0.0;
    double h = (ps[0]-ps[len-1])/len;
    result =  I_nu[0]*h*ps[0]/2.0;

    for (int idx = 1; idx < len-1; ++idx)
    {
        result += I_nu[idx]*h*ps[idx];
    }
    result += I_nu[len-1]*h*ps[len-1]/2.0;
    return result*8*M_PI*M_PI;
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

void integrate_source_functions(double* L_nu, const double* line_nu, const double* taus, const double* att_S_ul, const double* I_BB, 
        const double* nus, const double* ps, const double* Rs, double R_ph, double inv_ct, const int64_t* lens)
{
    double* I_nu  = calloc(lens[PLEN], sizeof(double));
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
                nu_lims = find_nu_limits_for_crossing_and_p(nus[nu_idx], ps[p_idx], cr_idx, no_of_cr_shells, inv_ct, Rs, line_nu, lens[LINELEN]);
                printf("Inner loop from %d to %d; ",nu_lims.start,nu_lims.end);
                for (int k_idx = nu_lims.start; k_idx < nu_lims.end; ++k_idx)
                {
                    I_nu[p_idx] = I_nu[p_idx] * exp(-taus[k_idx * lens[SHELLEN] + get_sh_idx(cr_idx,no_of_cr_shells)]) 
                                    + att_S_ul[k_idx * lens[SHELLEN] + get_sh_idx(cr_idx,no_of_cr_shells)];
                }
            }
        }
        L_nu[nu_idx] = integrate_intensity(I_nu, ps, lens[PLEN]); 
    }
    printf("\n\n");
}
