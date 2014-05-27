
#ifndef TARDIS_CMONTECARLO_H
#define TARDIS_CMONTECARLO_H

#include <Python.h>
#include <numpy/arrayobject.h>

#define MISS_DISTANCE 1e99

typedef struct StorageModel
{
  npy_float64 *packet_nus;
  npy_float64 *packet_mus;
  npy_float64 *packet_energies;
  npy_float64 *output_nus;
  npy_float64 *output_energies;
  npy_int64 *last_line_interaction_in_id;
  npy_int64 *last_line_interaction_out_id;
  npy_int64 *last_line_interaction_shell_id;
  npy_int64 *last_interaction_type;
  npy_int64 no_of_packets;
  npy_int64 no_of_shells;
  npy_float64 *r_inner;
  npy_float64 *r_outer;
  npy_float64 *v_inner;
  npy_float64 time_explosion;
  npy_float64 inverse_time_explosion;
  npy_float64 *electron_densities;
  npy_float64 *inverse_electron_densities;
  npy_float64 *line_list_nu;
  npy_float64 *line_lists_tau_sobolevs;
  npy_int64 line_lists_tau_sobolevs_nd;
  npy_float64 *line_lists_j_blues;
  npy_int64 line_lists_j_blues_nd;
  npy_int64 no_of_lines;
  npy_int64 line_interaction_id;
  npy_float64 *transition_probabilities;
  npy_int64 transition_probabilities_nd;
  npy_int64 *line2macro_level_upper;
  npy_int64 *macro_block_references;
  npy_int64 *transition_type;
  npy_int64 *destination_level_id;
  npy_int64 *transition_line_id;
  npy_float64 *js;
  npy_float64 *nubars;
  npy_float64 spectrum_start_nu;
  npy_float64 spectrum_delta_nu;
  npy_float64 spectrum_end_nu;
  npy_float64 *spectrum_virt_nu;
  npy_float64 sigma_thomson;
  npy_float64 inverse_sigma_thomson;
  npy_float64 inner_boundary_albedo;
  npy_int64 reflective_inner_boundary;
  npy_int64 current_packet_id;
} storage_model_t;

/** Look for a place to insert a value in an inversely sorted float array.
 *
 * @param x an inversely (largest to lowest) sorted float array
 * @param x_insert a value to insert
 * @param imin lower bound
 * @param imax upper bound
 *
 * @return index of the next boundary to the left
 */
npy_int64 binary_search(npy_float64 *x, npy_float64 x_insert, npy_int64 imin, npy_int64 imax);

/** Insert a value in to an array of line frequencies
 *
 * @param nu array of line frequencies
 * @param nu_insert value of nu key
 * @param number_of_lines number of lines in the line list
 *
 * @return index of the next line ot the red. If the key value is redder than the reddest line returns number_of_lines.
 */
npy_int64 line_search(npy_float64 *nu, npy_float64 nu_insert, npy_int64 number_of_lines);

/** Calculate the distance to the outer boundary.
 *
 * @return distance to the outer boundary
 */
inline npy_float64 compute_distance2outer(npy_float64 r, npy_float64 mu, npy_float64 r_outer);

/** Calculate the distance to the inner boundary.
 *
 * @return distance to the inner boundary
 */
inline npy_float64 compute_distance2inner(npy_float64 r, npy_float64 mu, npy_float64 r_inner);

#endif // TARDIS_CMONTECARLO_H
