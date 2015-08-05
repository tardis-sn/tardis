#ifndef TARDIS_STORAGE_H
#define TARDIS_STORAGE_H

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifdef __clang__
#define INLINE extern inline
#else
#define INLINE inline
#endif

typedef struct StorageModel
{
  double *packet_nus;
  double *packet_mus;
  double *packet_energies;
  double *output_nus;
  double *output_energies;
  int64_t *last_line_interaction_in_id;
  int64_t *last_line_interaction_out_id;
  int64_t *last_line_interaction_shell_id;
  int64_t *last_interaction_type;
  int64_t no_of_packets;
  int64_t no_of_shells;
  double *r_inner;
  double *r_outer;
  double *v_inner;
  double time_explosion;
  double inverse_time_explosion;
  double *electron_densities;
  double *inverse_electron_densities;
  double *line_list_nu;
  double *continuum_list_nu;
  double *line_lists_tau_sobolevs;
  int64_t line_lists_tau_sobolevs_nd;
  double *line_lists_j_blues;
  int64_t line_lists_j_blues_nd;
  int64_t no_of_lines;
  int64_t no_of_edges;
  int64_t line_interaction_id;
  double *transition_probabilities;
  int64_t transition_probabilities_nd;
  int64_t *line2macro_level_upper;
  int64_t *macro_block_references;
  int64_t *transition_type;
  int64_t *destination_level_id;
  int64_t *transition_line_id;
  double *transition_probabilities_continuum;
  int64_t transition_probabilities_nd_continuum;
  int64_t *cont_edge2macro_continuum;
  int64_t *macro_block_references_continuum;
  int64_t *transition_type_continuum;
  int64_t *destination_level_id_continuum;
  int64_t *transition_continuum_id;
  double *js;
  double *nubars;
  double spectrum_start_nu;
  double spectrum_delta_nu;
  double spectrum_end_nu;
  double spectrum_virt_start_nu;
  double spectrum_virt_end_nu;
  double *spectrum_virt_nu;
  double sigma_thomson;
  double inverse_sigma_thomson;
  double inner_boundary_albedo;
  int64_t reflective_inner_boundary;
  int64_t current_packet_id;
  double *chi_bf_tmp_partial;
  double *t_electrons;
  double *l_pop;
  double *l_pop_r;
  int64_t *ion_charge;
  double *ion_population;
  int64_t no_of_ions;
  ContinuumProcessesStatus cont_status;
  FreeFreeStatus ff_status;
} storage_model_t;

#endif // TARDIS_STORAGE_H
