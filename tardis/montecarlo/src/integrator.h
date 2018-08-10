#ifndef TARDIS_INTEGRATOR_H
#define TARDIS_INTEGRATOR_H

#include "storage.h"

double
integrate_intensity(const double* I_nu, const double h, int N);

int64_t
populate_z(const storage_model_t *storage, const double p, double *oz, int64_t *oshell_id);

void
calculate_p_values(double R_max, int64_t N, double *opp);

double
*_formal_integral(
                  const storage_model_t *storage, double T, double *nu, int64_t nu_size, double *att_S_ul, double *Jred_lu, double *Jblue_lu, int N);

#endif
