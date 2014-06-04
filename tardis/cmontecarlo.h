
#ifndef TARDIS_CMONTECARLO_H
#define TARDIS_CMONTECARLO_H

#include <Python.h>
#include <numpy/arrayobject.h>
#include "randomkit.h"

#define MISS_DISTANCE 1e99
#define C 29979245800.0
#define INVERSE_C 3.33564095198152e-11

rk_state mt_state;

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
 * @param r distance from the center to the packet
 * @param mu cosine of the angle the packet is moving at
 * @param r_outer distance from the center to the outer boundary
 *
 * @return distance to the outer boundary
 */
npy_float64 compute_distance2outer(npy_float64 r, npy_float64 mu, npy_float64 r_outer);

/** Calculate the distance to the inner boundary.
 *
 * @param r distance from the center to the packet
 * @param mu cosine of the angle the packet is moving at
 * @param r_inner distance from the center to the inner boundary
 *
 * @return distance to the inner boundary
 */
npy_float64 compute_distance2inner(npy_float64 r, npy_float64 mu, npy_float64 r_inner);

/** Calculate the distance the packet has to travel until it redshifts to the first spectral line.
 *
 * @param r distance from the center to the packet
 * @param mu angle at which the packet is moving
 * @param nu frequency of the packet
 * @param nu_line frequency of the next line the packet will encounter
 * @param t_exp time since the explosion
 * @param last_line 
 * @param next_line
 * @param cur_zone_id id of the current shell
 *
 * @return distance to the next spectral line
 */
npy_float64 compute_distance2line(npy_float64 r, npy_float64 mu, npy_float64 nu, npy_float64 nu_line, npy_float64 t_exp, npy_float64 inverse_t_exp, npy_float64 last_line, npy_float64 next_line, npy_int64 cur_zone_id);

/** Calculate the distance to the Thomson scatter event.
 * @param r distance from the center to the packet
 * @param mu cosine of the angle the packet is moving at
 * @param tau_event
 * @param inverse_ne
 *
 * @return distance to the Thomson scatter event in centimeters
 */
npy_float64 compute_distance2electron(npy_float64 r, npy_float64 mu, npy_float64 tau_event, npy_float64 inverse_ne);

inline npy_int64 macro_atom(npy_int64 activate_level, npy_float64 *p_transition, npy_int64 p_transition_nd, npy_int64 *type_transition, npy_int64 *target_level_id, npy_int64 *target_line_id, npy_int64 *unroll_reference, npy_int64 cur_zone_id);

inline npy_float64 move_packet(npy_float64 *r, npy_float64 *mu, npy_float64 nu, npy_float64 energy, npy_float64 distance, npy_float64 *js, npy_float64 *nubars, npy_float64 inverse_t_exp, npy_int64 cur_zone_id, npy_int64 virtual_packet);

inline void increment_j_blue_estimator(npy_int64 *current_line_id, npy_float64 *current_nu, npy_float64 *current_energy, npy_float64 *mu, npy_float64 *r, npy_float64 d_line, npy_int64 j_blue_idx, npy_float64 inverse_time_explosion, npy_float64 *line_lists_j_blues);

npy_int64 montecarlo_one_packet(storage_model_t *storage, npy_float64 *current_nu, npy_float64 *current_energy, npy_float64 *current_mu, npy_int64 *current_shell_id, npy_float64 *current_r, npy_int64 *current_line_id, npy_int64 *last_line, npy_int64 *close_line, npy_int64 *recently_crossed_boundary, npy_int64 virtual_packet_flag, npy_int64 virtual_mode);

npy_int64 montecarlo_one_packet_loop(storage_model_t *storage, npy_float64 *current_nu, npy_float64 *current_energy, npy_float64 *current_mu, npy_int64 *current_shell_id, npy_float64 *current_r, npy_int64 *current_line_id, npy_int64 *last_line, npy_int64 *close_line, npy_int64 *recently_crossed_boundary, npy_int64 virtual_packet_flag, npy_int64 virtual_packet);

#endif // TARDIS_CMONTECARLO_H
