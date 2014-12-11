#ifndef TARDIS_CMONTECARLO_H
#define TARDIS_CMONTECARLO_H

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "randomkit/randomkit.h"

#ifdef __clang__
#define INLINE extern inline
#else
#define INLINE inline
#endif

#define MISS_DISTANCE 1e99
#define C 29979245800.0
#define INVERSE_C 3.33564095198152e-11

typedef enum
{
  TARDIS_ERROR_OK = 0,
  TARDIS_ERROR_BOUNDS_ERROR = 1,
  TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE = 2
} tardis_error_t;

typedef enum
{
  TARDIS_PACKET_STATUS_IN_PROCESS = 0,
  TARDIS_PACKET_STATUS_EMITTED = 1,
  TARDIS_PACKET_STATUS_REABSORBED = 2
} rpacket_status_t;

/**
 * @brief A photon packet.
 */
typedef struct RPacket
{
  double nu; /**< Frequency of the packet in Hz. */
  double mu; /**< Cosine of the angle of the packet. */
  double energy; /**< Energy of the packet in erg. */
  double r; /**< Distance from center in cm. */
  double tau_event;
  double nu_line;
  int64_t current_shell_id; /**< ID of the current shell. */
  int64_t next_line_id;	/**< The index of the next line that the packet will encounter. */
  /**
   * @brief The packet has a nu red-ward of the last line.
   * It will not encounter any lines anymore.
   */
  int64_t last_line;
  /** 
   * @brief The packet just encountered a line that is very close to the next line.
   * The next iteration will automatically make an interaction with the next line 
   * (avoiding numerical problems).
   */
  int64_t close_line;
  /** 
   * @brief The packet has recently crossed the boundary and is now sitting on the boundary. 
   * To avoid numerical errors, make sure that d_inner is not calculated. The value is -1
   * if the packed moved inwards, 1 if the packet moved outwards and 0 otherwise.
   */
  int64_t recently_crossed_boundary;
  /**
   * @brief packet is a virtual packet and will ignore any d_line or d_electron checks.
   * It now whenever a d_line is calculated only adds the tau_line to an 
   * internal float.
   */
  int64_t virtual_packet_flag;
  int64_t virtual_packet;
  double d_line; /**< Distance to electron event. */
  double d_electron; /**< Distance to line event. */
  double d_boundary; /**< Distance to shell boundary. */
  int64_t next_shell_id; /**< ID of the next shell packet visits. */
  rpacket_status_t status; /**< Packet status (in process, emitted or reabsorbed). */
} rpacket_t;

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
  double *line_lists_tau_sobolevs;
  int64_t line_lists_tau_sobolevs_nd;
  double *line_lists_j_blues;
  int64_t line_lists_j_blues_nd;
  int64_t no_of_lines;
  int64_t line_interaction_id;
  double *transition_probabilities;
  int64_t transition_probabilities_nd;
  int64_t *line2macro_level_upper;
  int64_t *macro_block_references;
  int64_t *transition_type;
  int64_t *destination_level_id;
  int64_t *transition_line_id;
  double *js;
  double *nubars;
  double spectrum_start_nu;
  double spectrum_delta_nu;
  double spectrum_end_nu;
  double *spectrum_virt_nu;
  double sigma_thomson;
  double inverse_sigma_thomson;
  double inner_boundary_albedo;
  int64_t reflective_inner_boundary;
  int64_t current_packet_id;
} storage_model_t;

typedef void (*montecarlo_event_handler_t) (rpacket_t * packet,
					    storage_model_t * storage,
					    double distance);

/** Look for a place to insert a value in an inversely sorted float array.
 *
 * @param x an inversely (largest to lowest) sorted float array
 * @param x_insert a value to insert
 * @param imin lower bound
 * @param imax upper bound
 *
 * @return index of the next boundary to the left
 */
inline tardis_error_t reverse_binary_search (double *x, double x_insert,
					     int64_t imin, int64_t imax,
					     int64_t * result);

/** Look for a place to insert a value in a sorted float array.
 *
 * @param x a (lowest to largest) sorted float array
 * @param x_insert a value to insert
 * @param imin lower bound
 * @param imax upper bound
 *
 * @return index of the next boundary to the left
 */
inline tardis_error_t binary_search (double *x, double x_insert, int64_t imin,
				     int64_t imax, int64_t * result);

/** Insert a value in to an array of line frequencies
 *
 * @param nu array of line frequencies
 * @param nu_insert value of nu key
 * @param number_of_lines number of lines in the line list
 *
 * @return index of the next line ot the red. If the key value is redder than the reddest line returns number_of_lines.
 */
inline tardis_error_t line_search (double *nu, double nu_insert,
				   int64_t number_of_lines, int64_t * result);

/** Calculate the distance to shell boundary.
 *
 * @param packet rpacket structure with packet information
 * @param storage storage model data
 *
 * @return distance to shell boundary
 */
inline double compute_distance2boundary (rpacket_t * packet,
					 storage_model_t * storage);

/** Calculate the distance the packet has to travel until it redshifts to the first spectral line.
 *
 * @param packet rpacket structure with packet information
 * @param storage storage model data
 *
 * @return distance to the next spectral line
 */
extern inline tardis_error_t compute_distance2line (rpacket_t * packet,
						    storage_model_t * storage,
						    double *result);

/** Calculate the distance to the Thomson scatter event.
 *
 * @param packet rpacket structure with packet information
 * @param storage storage model data
 *
 * @return distance to the Thomson scatter event in centimeters
 */
inline double compute_distance2electron (rpacket_t * packet,
					 storage_model_t * storage);

inline int64_t macro_atom (rpacket_t * packet, storage_model_t * storage);

inline double move_packet (rpacket_t * packet, storage_model_t * storage,
			   double distance);

inline void increment_j_blue_estimator (rpacket_t * packet,
					storage_model_t * storage,
					double d_line, int64_t j_blue_idx);

int64_t montecarlo_one_packet (storage_model_t * storage, rpacket_t * packet,
			       int64_t virtual_mode);

int64_t montecarlo_one_packet_loop (storage_model_t * storage,
				    rpacket_t * packet,
				    int64_t virtual_packet);

inline double rpacket_get_nu (rpacket_t * packet);

inline void rpacket_set_nu (rpacket_t * packet, double nu);

inline double rpacket_get_mu (rpacket_t * packet);

inline void rpacket_set_mu (rpacket_t * packet, double mu);

inline double rpacket_get_energy (rpacket_t * packet);

inline void rpacket_set_energy (rpacket_t * packet, double energy);

inline double rpacket_get_r (rpacket_t * packet);

inline void rpacket_set_r (rpacket_t * packet, double r);

inline double rpacket_get_tau_event (rpacket_t * packet);

inline void rpacket_set_tau_event (rpacket_t * packet, double tau_event);

inline double rpacket_get_nu_line (rpacket_t * packet);

inline void rpacket_set_nu_line (rpacket_t * packet, double nu_line);

inline unsigned int rpacket_get_current_shell_id (rpacket_t * packet);

inline void rpacket_set_current_shell_id (rpacket_t * packet,
					  unsigned int current_shell_id);

inline unsigned int rpacket_get_next_line_id (rpacket_t * packet);

inline void rpacket_set_next_line_id (rpacket_t * packet,
				      unsigned int next_line_id);

inline bool rpacket_get_last_line (rpacket_t * packet);

inline void rpacket_set_last_line (rpacket_t * packet, bool last_line);

inline bool rpacket_get_close_line (rpacket_t * packet);

inline void rpacket_set_close_line (rpacket_t * packet, bool close_line);

inline int rpacket_get_recently_crossed_boundary (rpacket_t * packet);

inline void rpacket_set_recently_crossed_boundary (rpacket_t * packet,
						   int
						   recently_crossed_boundary);

inline int rpacket_get_virtual_packet_flag (rpacket_t * packet);

inline void rpacket_set_virtual_packet_flag (rpacket_t * packet,
					     int virtual_packet_flag);

inline int rpacket_get_virtual_packet (rpacket_t * packet);

inline void rpacket_set_virtual_packet (rpacket_t * packet,
					int virtual_packet);

inline double rpacket_get_d_boundary (rpacket_t * packet);

inline void rpacket_set_d_boundary (rpacket_t * packet, double d_boundary);

inline double rpacket_get_d_electron (rpacket_t * packet);

inline void rpacket_set_d_electron (rpacket_t * packet, double d_electron);

inline double rpacket_get_d_line (rpacket_t * packet);

inline void rpacket_set_d_line (rpacket_t * packet, double d_line);

inline int rpacket_get_next_shell_id (rpacket_t * packet);

inline void rpacket_set_next_shell_id (rpacket_t * packet, int next_shell_id);

inline rpacket_status_t rpacket_get_status (rpacket_t * packet);

inline void rpacket_set_status (rpacket_t * packet, rpacket_status_t status);

inline void rpacket_reset_tau_event (rpacket_t * packet);

tardis_error_t rpacket_init (rpacket_t * packet, storage_model_t * storage,
			     int packet_index, int virtual_packet_flag);

void initialize_random_kit (unsigned long seed);

#endif // TARDIS_CMONTECARLO_H
