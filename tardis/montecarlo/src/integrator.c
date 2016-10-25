#include <stdbool.h>
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

indexpair_t find_nu_limits_for_crossing_and_p(double nu, double p, int cr_idx, int no_of_cr_shells, double inv_ct, const double* Rs, const double* line_nu, int len, bool debug)
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
        nu_blu = nu * (1 - get_cr_sign(cr_idx,no_of_cr_shells)*z_blu*inv_ct);
        nu_red = nu * (1 - get_cr_sign(cr_idx+1,no_of_cr_shells)*z_red*inv_ct);
    }
    else 
    {
        if (p > Rs[1]) { 
            // Not intersecting photosphere
            z_cr = sqrt( Rs[0]*Rs[0] - p*p );
            nu_blu = nu * (1 + z_cr*inv_ct);
            nu_red = nu * (1 - z_cr*inv_ct); }
        else { 
            // Indexing follows from ordering and geometry
            nu_blu = nu * (1 - inv_ct*sqrt(Rs[1]*Rs[1]-p*p));
            nu_red = nu * (1 - inv_ct*sqrt(Rs[0]*Rs[0]-p*p));
        }

    }
    if ((nu_blu <= line_nu[len-1])|| (nu_red >= line_nu[0]) ) // Is nus in range of linelist? 
    {
        if(debug){
        printf("find_limits: no_cr %d ", no_of_cr_shells );
        printf(", nu_blu %e, nu_red %e, nu %e \n",nu_blu,nu_red, nu);}
        pair.start = -1; //Signal no line found
    }
    else 
    {
        if(debug){
        printf("find_limits: no_cr %d ", no_of_cr_shells );
        printf(", nu_blu %e, nu_red %e, nu %e ",nu_blu,nu_red, nu);
        printf("Nu in list\n");}
        pair = nu_idx_from_nu_pair(nu_blu, nu_red, line_nu, len);
    }

    return pair;
}
indexpair_t nu_idx_from_nu_pair(double nu_blu, double nu_red, const double* line_nu, int len)
{
    indexpair_t pair = { .start = 0, .end = len-1};
    bool found_blu = false;

    for (int idx = 0; idx < len; ++idx)
    {
        if ((line_nu[idx] <= nu_blu) & (!found_blu)){
            pair.start = idx;
            found_blu = true;}
        if (line_nu[idx] <= nu_red){
            pair.end = idx;
            break;}
    }

    return pair;
}

double intensity_black_body(double nu, double T){
    double k_B_cgs = 1.3806488000e-16;
    double h_cgs   = 6.6260695700e-27;
    double c_cgs   = 2.9979245800e+10;
    double coefficient, beta_rad;
    beta_rad = 1 / (k_B_cgs * T);
    coefficient = 2 * h_cgs / (c_cgs * c_cgs);
    return coefficient * nu*nu*nu / (exp(h_cgs * nu * beta_rad) -1 );
}    

double get_r(int cr_idx, int no_of_cr_shells, const double* Rs)
{
    int r_idx;
    if (cr_idx < no_of_cr_shells){
        r_idx = cr_idx;}
    else if (cr_idx < 2*no_of_cr_shells){
        r_idx = no_of_cr_shells-1 - cr_idx % no_of_cr_shells;}

    return Rs[r_idx];
}

int get_cr_sign(int cr_idx, int no_of_cr_shells)
{
    if (cr_idx < no_of_cr_shells){
        return -1;}
    else if (cr_idx < 2*no_of_cr_shells){
        return 1;}
}

int get_cr_start(int no_of_shells, double p, double R_ph)
{
    if (p >= R_ph) {
        return 0;}
    else if (p < R_ph) {
        return no_of_shells+1;
    }
}

int get_sh_idx(int cr_idx, int no_of_cr_shells)
{
    if (cr_idx < no_of_cr_shells){
       return cr_idx;}
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
double sum_lines(indexpair_t nu_lims, double I_nu, const double* taus, const double* att_S_ul, int sh_idx, int len)
{
    double result = I_nu;
    for (int k_idx = nu_lims.start; k_idx <= nu_lims.end; ++k_idx)
    {
        result = result * exp(-taus[k_idx*len +sh_idx]) + GET_IJ(att_S_ul,sh_idx,k_idx,len);
    }
    return result;
}
double integrate_intensity(const double* I_nu, const double* ps, int len)
{
    double result = 0.0;
    double h = (ps[0]-ps[len-1])/(len-1);
    result =  I_nu[0]*h*ps[0]/2.0;

    for (int idx = 1; idx < len-1; ++idx)
    {
        result += I_nu[idx]*h*ps[idx];
    }
    result += I_nu[len-1]*h*ps[len-1]/2.0;
    return result*8*M_PI*M_PI;
}
void integrate_source_functions(double* L_nu, const double* line_nu, const double* taus, const double* att_S_ul, const double* nus,
        const double* ps, const double* Rs, double T, double R_ph, double inv_ct, const int64_t* lens)
{
    double* I_nu  = calloc(lens[PLEN], sizeof(double));
    int no_of_cr_shells, sh_idx, cr_start, cr_end;
    int kludge = 0;
    indexpair_t nu_lims;
    for (int nu_idx = 0; nu_idx < lens[NULEN]; ++nu_idx)
    {
        memset(I_nu,0.0, lens[PLEN] * sizeof(I_nu));
        for (int p_idx = 0; p_idx < lens[PLEN]; ++p_idx)
        {
            no_of_cr_shells = get_num_shell_cr(ps[p_idx],Rs,lens[SHELLEN]);
            if (ps[p_idx] > R_ph) 
            {
                kludge   = 0;
                cr_start = 0;
                cr_end   = 2*no_of_cr_shells;
            }
            else
            {
                I_nu[p_idx] = intensity_black_body(nus[nu_idx], T);
                cr_start = lens[SHELLEN]+1; // get_cr_start(lens[SHELLEN], ps[p_idx], R_ph);
                cr_end   = 2*no_of_cr_shells+1;
                kludge   = 1; // Needed to get proper sh_idx when photosphere is crossed
                
                find_nu_limits_for_crossing_and_p(nus[nu_idx], ps[p_idx], 2, no_of_cr_shells, inv_ct, Rs, line_nu, lens[LINELEN],true);
            }

            for (int cr_idx = cr_start; cr_idx < cr_end; ++cr_idx)
            {
                nu_lims = find_nu_limits_for_crossing_and_p(nus[nu_idx], ps[p_idx], cr_idx, no_of_cr_shells, inv_ct, Rs, line_nu, lens[LINELEN],false);
                sh_idx  = get_sh_idx(cr_idx,no_of_cr_shells+kludge);
                if (nu_lims.start > -1) { // Just sum lines if any line was found                    
                    I_nu[p_idx] = sum_lines(nu_lims, I_nu[p_idx], taus, att_S_ul, sh_idx, lens[SHELLEN]);}
            }
        }
        L_nu[nu_idx] = integrate_intensity(I_nu, ps, lens[PLEN]); 
    }
}


double test_nu_limits_for_crossing_and_p(double nu, double p, int cr_idx, int no_of_cr_shells, double inv_ct, const double* Rs, const double* line_nu, int len, int red)
{
    double blu_R, red_R, z_blu, z_red, z_cr, nu_blu, nu_red;
    double out;

    if (no_of_cr_shells > 1)
    {
        assert(Rs[0] > Rs[1]); // Decreasing order
        
        blu_R = get_r(cr_idx,no_of_cr_shells,Rs);
        red_R = get_r(cr_idx+1,no_of_cr_shells,Rs);
        z_blu = sqrt( blu_R*blu_R - p*p );
        z_red = sqrt( red_R*red_R - p*p );
        nu_blu = nu * (1 - get_cr_sign(cr_idx,no_of_cr_shells)*z_blu*inv_ct);
        nu_red = nu * (1 - get_cr_sign(cr_idx+1,no_of_cr_shells)*z_red*inv_ct);
    }
    else 
    {
        z_cr = sqrt( Rs[cr_idx]*Rs[cr_idx] - p*p );
        nu_blu = nu * (1 - z_cr*inv_ct);
        nu_red = nu * (1 + z_cr*inv_ct);
    }
    out  = nu_blu;
    if (red != 0){
    out  = nu_red;}
    return out;
}

double test_index_macro(const double* array, int i, int j, int jlen)
{
    return GET_IJ(array,i,j,jlen);
}
void debug_print_arg(double* arg,int len)
{
    for (int64_t i = 0; i < len; i++)
    {
        printf("%d: %f, ", i, arg[i]);
        if ( (i % 8 == 0 ) & i > 0){printf("\n");}
    }
}
void debug_print_2d_arg(double* arg,int len1, int len2)
{
    for (int64_t j = 0; j < len1; ++j)
    {
        for(int64_t i = 0; i < len2; ++i)
        {
        printf("[%d,%d:%d]: %.8f, ",i,j,i*len1+j,arg[i*len1+j]);
        }
        printf("\n\n");
    } 
}

