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
#include "cmontecarlo.h"
#include "cmontecarlo1.h"

double calc_f(storage_model_t * storage, rpacket_t * packet) // Calcuates value for the filling factor f from velocity
{
  double f_0 = storage->filling_factor; // Initial filling factor.
  double R = storage->density_ratio; // Ratio of cloud density to intercloud density
  double v_initial = 1.10000000e+09; // Initial velocity
  double k = storage->power_law_k; // Constant for power law calculation
  double v_in = storage->v_inner[rpacket_get_current_shell_id (packet)];

  double f = f_0*pow((v_in/v_initial),k);

  return f;
}

double cloud_tau(storage_model_t *storage, double tau_line, double f) // Calculates tau_sobolev for photons in cloud
{
  double R = storage->density_ratio;

  tau_line = tau_line*(1/(f+((1-f)*1/R))); // Cloud optical depth

  return tau_line;
}

double intercloud_tau(storage_model_t *storage, double tau_line, double f) // Calculates tau_sobolev for photons outside cloud
{
  double R = storage->density_ratio;

  tau_line = tau_line*(1/((f*R)+(1-f))); // Intercloud optical depth

  return tau_line;
}

double calc_tau(storage_model_t *storage, rpacket_t *packet, int64_t line2d_idx, rk_state *mt_state)
{
  double random_clump = rk_double(mt_state); // Generates a random number to determine if photon is in cloud.

  double f = calc_f(storage, packet); // Calculate volume filling factor as a function of ejecta velocity.

  double tau_line =
    storage->line_lists_tau_sobolevs[line2d_idx];

    if(random_clump<f)
      {
        tau_line = cloud_tau(storage, tau_line, f); // Cloud optical depth
      }
    else
      {
        tau_line = intercloud_tau(storage, tau_line, f); // Intercloud optical depth
      }

  return tau_line;
}

#endif // CLUMPING.H
