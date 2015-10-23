#ifndef clumping.h
#define clumping.h

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "randomkit/randomkit.h"
#include "status.h"
#include "storage.h"
#include "cmontecarlo1.h"

double calc_f(storage_model_t *storage, double v_in)
{
  double f_0 = storage->filling_factor; // True "clumping" factor.
  double R = storage->density_ratio; // Ratio of cloud density to intercloud density
  double v_initial = 1.10000000e+09;
  double k = storage->power_law_k;

  double f = f_0*pow((v_in/v_initial),k); // Calculate volume filling factor as a function of ejecta velocity.

  return f;
}

double cloud_tau(storage_model_t *storage, double tau_line, double f)
{
  double R = storage->density_ratio;

  tau_line = tau_line*(1/(f+((1-f)*1/R)));

  return tau_line
}

double intercloud_tau(storage_model_t *storage, double tau_line, double f)
{
  double R = storage->density_ratio;

  tau_line = tau_line*(1/((f*R)+(1-f))); // Intercloud optical depth

  return tau_line
}

#endif // CLUMPING.H
