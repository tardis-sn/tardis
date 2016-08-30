#ifndef TARDIS_INTEGRATOR_H
#define TARDIS_INTEGRATOR_H

#include "storage.h"

typedef struct IndexPair
{
    int start;
    int end;
} indexpair_t;

int get_cr_sign(int cr_idx, int no_of_shell_cr);
int get_cr_start(int no_of_shell_cr, double p, double R_ph);
int get_sh_idx(int cr_idx, int no_of_shell_cr);
int get_num_shell_cr(double p, double* Rs, int len);

double get_r(int cr_idx, int no_of_shell_cr, double* Rs);
double integrate_intensity(double* I_nu, double* ps, int len);

void integrate_source_functions(double* L_nu, double* line_nu, double* taus, double* att_S_ul, double* I_BB, double* nus, 
              double* ps, double* Rs, double R_ph, int64_t* lens);

indexpair_t find_nu_limits_for_crossing_and_p(double nu, double p,int cr_idx, int no_of_shell_cr, double* Rs, double inv_t);

void debug_print_arg(double* arg,int len);
void debug_print_2d_arg(double* arg,int len1, int len2);

#endif
