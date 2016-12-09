#ifndef TARDIS_INTEGRATOR_H
#define TARDIS_INTEGRATOR_H

#include "storage.h"

typedef struct IndexPair
{
    int start;
    int end;
} indexpair_t;

indexpair_t find_nu_limits_for_crossing_and_p(double nu, double p, int cr_idx, int no_of_cr_shells, double inv_t, const double* Rs, const double* line_nu, int len);

int get_cr_sign(int cr_idx, int no_of_cr_shells);
int get_cr_start(int no_of_cr_shells, double p, double R_ph);
int get_sh_idx(int cr_idx, int no_of_cr_shells);
int get_num_shell_cr(double p, const double* Rs, int len);

double get_r(int cr_idx, int no_of_cr_shells, const double* Rs);
double integrate_intensity(const double* I_nu, const double h, int N);

void debug_print_arg(double* arg,int len);
void debug_print_2d_arg(double* arg,int len1, int len2);
void integrate_source_functions(double* L_nu, const double* line_nu, const double* taus, const double* att_S_ul, const double* I_BB, 
        const double* nus, const double* ps, const double* Rs, double R_ph, const int64_t* lens);

void _formal_integral(
        storage_model_t *storage, double *I_BB, double *att_S_ul, int N, double *L);
/*
double recursive_update_luminosity(
    double nu_start, double nu_end, double *line_nu, double *tau, double *att_S_ul, double initial);
    */

#endif
