#ifndef TARDIS_INTEGRATOR_H
#define TARDIS_INTEGRATOR_H

#include "storage.h"

typedef struct IndexPair
{
    int start;
    int end;
} indexpair_t;

void integrate_source_functions(double* L_nu, double* line_nu, double* taus, double* att_S_ul, double* I_BB, double* nus, 
              double* ps, double* Rs, double R_ph, double ct, int* lens);

double integrate_intensity(double* I_nu,double* ps, int* len);

indexpair_t find_array_bounds_nonzero(double* array,int idx2, int len);
indexpair_t find_nu_limits_for_shell_and_p(double p,int shl_idx, double* Rs, int no_shells);

void debug_print_arg(double* arg,int len);

void debug_print_2d_arg(double* arg,int len1, int len2);

#endif
