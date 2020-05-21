#define _USE_MATH_DEFINES

#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "io.h"
#include "storage.h"
#include "integrator.h"
#include "cmontecarlo.h"

#include "omp_helper.h"

#define NULEN   0
#define LINELEN 1
#define PLEN    2
#define SHELLEN 3

#define C_INV 3.33564e-11
#define M_PI acos (-1)
#define KB_CGS 1.3806488e-16
#define H_CGS 6.62606957e-27

/**
 * Calculate the intensity of a black-body according to the following formula
 * .. math::
 * I(\\nu, T) = \\frac{2h\\nu^3}{c^2}\frac{1}{e^{h\\nu \\beta_\\textrm{rad}} - 1}
*/
double
intensity_black_body (double nu, double T)
{
  double beta_rad = 1 / (KB_CGS * T);
  double coefficient = 2 * H_CGS * C_INV * C_INV;
  return coefficient * nu * nu * nu / (exp(H_CGS * nu * beta_rad) - 1 );
}


/*! @brief Algorithm to integrate an array using the trapezoid integration rule
 *
*/
double
trapezoid_integration (const double* array, const double h, int N)
{
  double result = (array[0] + array[N-1])/2;
  for (int idx = 1; idx < N-1; ++idx)
    {
      result += array[idx];
    }
  return result * h;
}

/*! @brief Calculate distance to p line
 *
 *  Calculate half of the length of the p-line inside a shell
 *  of radius r in terms of unit length (c * t_exp).
 *  If shell and p-line do not intersect, return 0.
 *
 * @param r radius of the shell
 * @param p distance of the p-line to the center of the supernova
 * @param inv_t inverse time_explosio is needed to norm to unit-length
 * @return half the lenght inside the shell or zero
 */
static inline double
calculate_z(double r, double p, double inv_t)
{
  return (r > p) ? sqrt(r * r - p * p) * C_INV * inv_t : 0;
}


/*!
 * @brief Calculate p line intersections
 *
 * This function calculates the intersection points of the p-line with each shell
 *
 * @param storage (INPUT) A storage model containing the environment
 * @param p (INPUT) distance of the integration line to the center
 * @param oz (OUTPUT) will be set with z values. The array is truncated by the
 *                  value `1`.
 * @param oshell_id (OUTPUT) will be set with the corresponding shell_ids
 * @return number of shells intersected by the p-line
 */
int64_t
populate_z(const storage_model_t *storage, const double p, double *oz, int64_t *oshell_id)
{

  // Abbreviations
  double *r = storage->r_outer_i;
  const int64_t N = storage->no_of_shells_i;
  double inv_t = storage->inverse_time_explosion;
  double z = 0;

  int64_t i = 0, offset = N, i_low, i_up;

  if (p <= storage->r_inner_i[0])
    {
      // Intersect the photosphere
      for(i = 0; i < N; ++i)
        { // Loop from inside to outside
          oz[i] = 1 - calculate_z(r[i], p, inv_t);
          oshell_id[i] = i;
        }
      return N;
    }
  else
    {
      // No intersection with the photosphere
      // that means we intersect each shell twice
      for(i = 0; i < N; ++i)
        { // Loop from inside to outside
          z = calculate_z(r[i], p, inv_t);
          if (z == 0)
            continue;
          if (offset == N)
            {
              offset = i;
            }
          // Calculate the index in the resulting array
          i_low = N - i - 1;  // the far intersection with the shell
          i_up = N + i - 2 * offset; // the nearer intersection with the shell

          // Setting the arrays
          oz[i_low] = 1 + z;
          oshell_id[i_low] = i;
          oz[i_up] = 1 - z;
          oshell_id[i_up] = i;
        }
      return 2 * (N - offset);
    }
}


/*! @brief Calculate integration points
 *
 */
void
calculate_p_values(double R_max, int64_t N, double *opp)
{
  for(int i = 0; i<N; ++i)
    {
      // Trapezoid integration points
      opp[i] = R_max/(N - 1) * (i);
    }
}

/*! @brief Caculate a spectrum using the formal integral approach
 *
 */
double *
_formal_integral(
                 const storage_model_t *storage,
                 double iT,
                 double *inu, int64_t inu_size,
                 double *att_S_ul, double *Jred_lu, double *Jblue_lu, int N)
{

  // Initialize the output which is shared among threads
  double *L = calloc(inu_size, sizeof(double));

  // global read-only values
  int64_t size_line = storage->no_of_lines,
          size_shell = storage->no_of_shells_i,
          size_tau = size_line * size_shell,
          finished_nus = 0;

  double R_ph = storage->r_inner_i[0];
  double R_max = storage->r_outer_i[size_shell - 1];
  double pp[N];
  double *exp_tau = calloc(size_tau, sizeof(double));
//#pragma omp parallel firstprivate(L, exp_tau)
#pragma omp parallel 
    {

#pragma omp master
        {
          if (omp_get_num_threads() > 1) {
              fprintf(stderr, "Doing the formal integral\nRunning with OpenMP - %d threads\n", omp_get_num_threads());
          } else {
              fprintf(stderr, "Doing the formal integral\nRunning without OpenMP\n");
          }
          print_progress_fi(0, inu_size);
        }

      // Initializing all the thread-local variables
      int64_t offset = 0, i = 0,
              size_z = 0,
              idx_nu_start = 0,
              direction = 0,
              first = 0;

      double I_nu[N],
             //I_nu_b[N],
             //I_nu_r[N],
             z[2 * storage->no_of_shells_i],
             p = 0,
             nu_start,
             nu_end,
             nu,
             zstart,
             zend,
             escat_contrib,
             escat_op,
             Jkkp;
      int64_t shell_id[2 * storage->no_of_shells_i];

      double *pexp_tau, *patt_S_ul, *pline, *pJred_lu, *pJblue_lu;

      // Prepare exp_tau
#pragma omp for
      for (i = 0; i < size_tau; ++i) {
          exp_tau[i] = exp( -storage->line_lists_tau_sobolevs_i[i]);
      }
      calculate_p_values(R_max, N, pp);
      // Done with the initialization

      // Loop over wavelengths in spectrum
#pragma omp for
      for (int nu_idx = 0; nu_idx < inu_size ; ++nu_idx)
        {
          nu = inu[nu_idx];

          // Loop over discrete values along line
          for (int p_idx = 1; p_idx < N; ++p_idx)
            {
              escat_contrib = 0;
              p = pp[p_idx];

              // initialize z intersections for p values
              size_z = populate_z(storage, p, z, shell_id);

              // initialize I_nu
              if (p <= R_ph)
                I_nu[p_idx] = intensity_black_body(nu * z[0], iT);
              else
                I_nu[p_idx] = 0;

              // Find first contributing line
              nu_start = nu * z[0];
              nu_end = nu * z[1];
              line_search(
                  storage->line_list_nu,
                  nu_start,
                  size_line,
                  &idx_nu_start
              );
              offset = shell_id[0] * size_line;

              // start tracking accumulated e-scattering optical depth
              zstart = storage->time_explosion / C_INV * (1. - z[0]);

              // Initialize pointers
              pline = storage->line_list_nu + idx_nu_start;
              pexp_tau = exp_tau + offset + idx_nu_start;
              patt_S_ul = att_S_ul + offset + idx_nu_start;
              pJred_lu = Jred_lu + offset + idx_nu_start;
              pJblue_lu = Jblue_lu + offset + idx_nu_start;

              // flag for first contribution to integration on current p-ray
              first = 1;

              // TODO: Ugly loop
              // Loop over all intersections
              // TODO: replace by number of intersections and remove break
              for (i = 0; i < size_z - 1; ++i)
                {
                  escat_op = storage->electron_densities_i[shell_id[i]] * storage->sigma_thomson;
                  nu_end = nu * z[i+1];

                  // TODO: e-scattering: in principle we also have to check
                  // that dtau is <<1 (as assumed in Lucy 1999); if not, there
                  // is the chance that I_nu_b becomes negative
                  for (;pline < storage->line_list_nu + size_line;
                       // We have to increment all pointers simultaneously
                       ++pline,
                       ++pexp_tau,
                       ++patt_S_ul,
                       ++pJblue_lu)
                    {
                      if (*pline < nu_end)
                      {
                        // next resonance not in current shell
                        break;
                      }

                      // Calculate e-scattering optical depth to next resonance point
                      zend = storage->time_explosion / C_INV * (1. - *pline / nu);

                      if (first == 1){
                        // First contribution to integration
                        // NOTE: this treatment of I_nu_b (given by boundary
                        // conditions) is not in Lucy 1999; should be
                        // re-examined carefully 
                        escat_contrib += (zend - zstart) * escat_op * (*pJblue_lu - I_nu[p_idx]) ;
                        first = 0;
                      }
                      else{
                        // Account for e-scattering, c.f. Eqs 27, 28 in Lucy 1999
                        Jkkp = 0.5 * (*pJred_lu + *pJblue_lu);
                        escat_contrib += (zend - zstart) * escat_op * (Jkkp - I_nu[p_idx]) ;
                        // this introduces the necessary offset of one element between pJblue_lu and pJred_lu
                        pJred_lu += 1;
                      }
                      I_nu[p_idx] = I_nu[p_idx] + escat_contrib;

                      // Lucy 1999, Eq 26
                      I_nu[p_idx] = I_nu[p_idx] * (*pexp_tau) + *patt_S_ul;

                      // reset e-scattering opacity 
                      escat_contrib = 0;
                      zstart = zend;
                    }
                    // Calculate e-scattering optical depth to grid cell boundary
                    Jkkp = 0.5 * (*pJred_lu + *pJblue_lu);
                    zend = storage->time_explosion / C_INV * (1. - nu_end / nu);
                    escat_contrib += (zend - zstart) * escat_op * (Jkkp - I_nu[p_idx]);
                    zstart = zend;

                    if (i < size_z-1){
                      // advance pointers
                      direction = shell_id[i+1] - shell_id[i];
                      pexp_tau  += direction * size_line;
                      patt_S_ul += direction * size_line;
                      pJred_lu  += direction * size_line;
                      pJblue_lu += direction * size_line;
                    }
                }
              I_nu[p_idx] *= p;
            }
          // TODO: change integration to match the calculation of p values
          L[nu_idx] = 8 * M_PI * M_PI * trapezoid_integration(I_nu, R_max/N, N);
#pragma omp atomic update
          ++finished_nus;
          if (finished_nus%10 == 0){
              print_progress_fi(finished_nus, inu_size);
          }
        }
      // Free everything allocated on heap
      printf("\n");
    }
      free(exp_tau);
  return L;
}
