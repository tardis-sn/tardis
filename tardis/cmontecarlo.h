#ifndef TARDIS_CMONTECARLO_H
#define TARDIS_CMONTECARLO_H

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "randomkit.h"

#define MISS_DISTANCE 1e99
#define C 29979245800.0
#define INVERSE_C 3.33564095198152e-11
#define H 6.6260755e-27		// erg*s, converted to CGS units from the NIST Constant Index
#define KB 1.3806488e-16	//erg / K converted to CGS units from the NIST Constant Index
#define BF_DELTA_NU_RECOMPUTE 100 // Allowed deviation from the nu used for the bf chi until the stored chi is recomputed
#define TRUE 1

#define FALSE 0
/*
defines for the packet_status_t flags
[R,K,I,IN_PROCESS,STATUS_EMITTED,STATUS_REABSORBED,STATUS_DISABLED]
*/
#define R_PACKET 1 << 6
#define K_PACKET 1 << 5
#define I_PACKET 1 << 4
#define PACKET_IN_PROCESS 1 << 3
#define PACKET_STATUS_DISABLED 1 << 0
#define PACKET_STATUS_REABSORBED 1 << 1
#define PACKET_STATUS_EMITTED 1 << 2



typedef enum {
	TARDIS_ERROR_OK = 0,
	TARDIS_ERROR_BOUNDS_ERROR = 1,
	TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE = 2
} tardis_error_t;

typedef enum			// This is the new main status for packet types
{
	//[IN_PROCESS,STATUS_EMITTED,STATUS_REABSORBED,STATUS_DISABLED]
	TARDIS_PACKET_STATUS_IN_PROCESS = PACKET_IN_PROCESS,
	TARDIS_PACKET_STATUS_EMITTED = PACKET_STATUS_EMITTED,
	TARDIS_PACKET_STATUS_REABSORBED = PACKET_STATUS_REABSORBED,
	TARDIS_PACKET_STATUS_DISABLED = PACKET_STATUS_DISABLED
} rpacket_status_t;

typedef enum {
	// > 0 IN_PROCESS; > 10 EMITTED ; > 20 REABSORBED ; > 30 DISABLED
	// [R,K,I,IN_PROCESS,STATUS_EMITTED,STATUS_REABSORBED,STATUS_DISABLED]
	TARDIS_R_PACKET_IN_PROCESS = R_PACKET | PACKET_IN_PROCESS,
	TARDIS_R_PACKET_STATUS_EMITTED = R_PACKET | PACKET_STATUS_EMITTED,
	TARDIS_R_PACKET_STATUS_REABSORBED = R_PACKET | PACKET_STATUS_REABSORBED,
	TARDIS_R_PACKET_STATUS_DISABLED = R_PACKET | PACKET_STATUS_DISABLED,
	TARDIS_K_PACKET_IN_PROCESS = K_PACKET | PACKET_IN_PROCESS,
	TARDIS_K_PACKET_STATUS_DISABLED = K_PACKET | PACKET_STATUS_DISABLED,
	TARDIS_I_PACKET_IN_PROCESS = I_PACKET | PACKET_IN_PROCESS,
	TARDIS_I_PACKET_STATUS_DISABLED = I_PACKET | PACKET_STATUS_DISABLED
} packet_status_t;

/**
 * @brief A photon packet.
 */
typedef struct RPacket {
	double nu;
	     /**< Frequency of the packet in Hz. */
	double mu;
	     /**< Cosine of the angle of the packet. */
	double energy;
		 /**< Energy of the packet in erg. */
	double comov_nu;
		   /**<Frequency of the packet in the comoving frame. */
	double comov_energy;
		       /**< Energy of the packet in the comoving frame. */
	double r;
	    /**< Distance from center in cm. */
	double tau_event;
		    /**< Optical depth to next event. */
	double nu_line;
		  /**< frequency of the last line. */
	int64_t current_shell_id;
			    /**< ID of the current shell. */
	int64_t next_line_id;
			/**< The index of the next line that the packet will encounter. */
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
	double d_inner;
		  /**< Distance to the inner shell boundary. */
	double d_outer;
		  /**< Distance to the outer shell boundary. */
	double d_line;
		 /**< Distance to the next possible line event. */
	double d_electron;
		     /**< Distance to the next electron scatter event. */
	int64_t moved;
	double d_boundary;
		     /**< Distance to shell boundary. */
	int64_t next_shell_id;
			 /**< ID of the next shell packet visits. */
	rpacket_status_t status;
			   /**< Packet status (in process, emitted or reabsorbed). */
	packet_status_t packet_status;
	double chi_bf;
	double chi_th;
	double chi_ff;
	double chi_cont;
	double d_bf;
	double d_th;
	double d_ff;
	double d_cont;
	double last_bf_edge;
	double *chi_bf_partial;
	int64_t chi_bf_tmp_partial_last_shell_id;
	double chi_bf_tmp_partial_last_nu;
} rpacket_t;

typedef struct StorageModel {
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

	int64_t *chi_bf_index_to_level;
	int64_t chi_bf_index_to_level_nrow;
	int64_t chi_bf_index_to_level_ncolum;

	double *bf_level_population;
	int64_t bf_level_population_nrow;
	int64_t bf_level_population_ncolum;

	double *bf_lpopulation_ratio;
	int64_t bf_lpopulation_ratio_nrow;
	int64_t bf_lpopulation_ratio_ncolum;

	double *bf_lpopulation_ratio_nlte_lte;
	int64_t bf_lpopulation_ratio_nlte_lte_nrow;
	int64_t bf_lpopulation_ratio_nlte_lte_ncolum;

	double *bf_cross_sections;
	double *bound_free_th_frequency;

	double *t_electrons;
	double *chi_bf_tmp_partial;


    double *Cr_fb_kji_all;
    int64_t Cr_fb_kji_all_nrow;
    int64_t Cr_fb_kji_all_ncolum;

    double *Cr_fb_kji_cumsum_all;
    int64_t Cr_fb_kji_cumsum_all_nrow;
    int64_t Cr_fb_kji_cumsum_all_ncolum;

    int64_t *Cr_fb_kji_index;


    double *Cr_fb_kj_all;
    int64_t Cr_fb_kj_all_nrow;
    int64_t Cr_fb_kj_all_ncolum;

    double *Cr_fb_kj_cumsum_all;
    int64_t Cr_fb_kj_cumsum_all_nrow;
    int64_t Cr_fb_kj_cumsum_all_ncolum;

    int64_t *Cr_fb_kj_index;


    double *Cr_bb_kij_all;
    int64_t Cr_bb_kij_all_nrow;
    int64_t Cr_bb_kij_all_ncolum;

    double *Cr_bb_kij_cumsum_all;
    int64_t Cr_bb_kij_cumsum_all_nrow;
    int64_t Cr_bb_kij_cumsum_all_ncolum;

    int64_t *Cr_bb_kij_index;
	//double kB;

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
inline tardis_error_t reverse_binary_search(double *x, double x_insert, int64_t imin,
				    int64_t imax, int64_t * result);



/** Look for a place to insert a value in a sorted float array.
 *
 * @param x a (lowest to largest) sorted float array
 * @param x_insert a value to insert
 * @param imin lower bound
 * @param imax upper bound
 *
 * @return index of the next boundary to the left
 */
inline tardis_error_t binary_search(double *x, double x_insert, int64_t imin,
				    int64_t imax, int64_t * result);

/** Insert a value in to an array of line frequencies
 *
 * @param nu array of line frequencies
 * @param nu_insert value of nu key
 * @param number_of_lines number of lines in the line list
 *
 * @return index of the next line ot the red. If the key value is redder than the reddest line returns number_of_lines.
 */
inline tardis_error_t line_search(double *nu, double nu_insert,
				  int64_t number_of_lines, int64_t * result);

/** Calculate the distance to shell boundary.
 *
 * @param packet rpacket structure with packet information
 * @param storage storage model data
 *
 * @return distance to shell boundary
 */
inline double compute_distance2boundary(rpacket_t * packet,
					storage_model_t * storage);

/** Calculate the distance the packet has to travel until it redshifts to the first spectral line.
 *
 * @param packet rpacket structure with packet information
 * @param storage storage model data
 *
 * @return distance to the next spectral line
 */
extern inline tardis_error_t compute_distance2line(rpacket_t * packet,
						   storage_model_t * storage,
						   double *result);

/** Calculate the distance to the Thomson scatter event.
 *
 * @param packet rpacket structure with packet information
 * @param storage storage model data
 *
 * @return distance to the Thomson scatter event in centimeters
 */
inline double compute_distance2electron(rpacket_t * packet,
					storage_model_t * storage);

inline int64_t macro_atom(rpacket_t * packet, storage_model_t * storage);

inline double move_packet(rpacket_t * packet, storage_model_t * storage,
			  double distance);

inline void increment_j_blue_estimator(rpacket_t * packet,
				       storage_model_t * storage, double d_line,
				       int64_t j_blue_idx);

int64_t montecarlo_one_packet(storage_model_t * storage, rpacket_t * packet,
			      int64_t virtual_mode);

int64_t montecarlo_one_packet_loop(storage_model_t * storage,
				   rpacket_t * packet, int64_t virtual_packet);

inline double rpacket_get_nu(rpacket_t * packet);

inline void rpacket_set_nu(rpacket_t * packet, double nu);

inline double rpacket_get_mu(rpacket_t * packet);

inline void rpacket_set_mu(rpacket_t * packet, double mu);

inline double rpacket_get_energy(rpacket_t * packet);

inline void rpacket_set_energy(rpacket_t * packet, double energy);

inline double rpacket_get_r(rpacket_t * packet);

inline void rpacket_set_r(rpacket_t * packet, double r);

inline double rpacket_get_tau_event(rpacket_t * packet);

inline void rpacket_set_tau_event(rpacket_t * packet, double tau_event);

inline double rpacket_get_nu_line(rpacket_t * packet);

inline void rpacket_set_nu_line(rpacket_t * packet, double nu_line);

inline unsigned int rpacket_get_current_shell_id(rpacket_t * packet);

inline void rpacket_set_current_shell_id(rpacket_t * packet,
					 unsigned int current_shell_id);

inline unsigned int rpacket_get_next_line_id(rpacket_t * packet);

inline void rpacket_set_next_line_id(rpacket_t * packet,
				     unsigned int next_line_id);

inline bool rpacket_get_last_line(rpacket_t * packet);

inline void rpacket_set_last_line(rpacket_t * packet, bool last_line);

inline bool rpacket_get_close_line(rpacket_t * packet);

inline void rpacket_set_close_line(rpacket_t * packet, bool close_line);

inline int rpacket_get_recently_crossed_boundary(rpacket_t * packet);

inline void rpacket_set_recently_crossed_boundary(rpacket_t * packet,
						  int
						  recently_crossed_boundary);

inline int rpacket_get_virtual_packet_flag(rpacket_t * packet);

inline void rpacket_set_virtual_packet_flag(rpacket_t * packet,
					    int virtual_packet_flag);

inline int rpacket_get_virtual_packet(rpacket_t * packet);

inline void rpacket_set_virtual_packet(rpacket_t * packet, int virtual_packet);

inline void rpacket_set_d_boundary(rpacket_t * packet, double d_boundary);

inline double rpacket_get_d_electron(rpacket_t * packet);

inline void rpacket_set_d_electron(rpacket_t * packet, double d_electron);

inline double rpacket_get_d_line(rpacket_t * packet);

inline void rpacket_set_d_line(rpacket_t * packet, double d_line);

inline int rpacket_get_next_shell_id(rpacket_t * packet);

inline void rpacket_set_next_shell_id(rpacket_t * packet, int next_shell_id);

inline double rpacket_get_d_boundary(rpacket_t * packet);

inline void rpacket_set_d_boundary(rpacket_t * packet, double d_boundary);

inline double rpacket_get_d_continuum(rpacket_t * packet);

inline double rpacket_get_d_electron(rpacket_t * packet);

inline double rpacket_get_d_freefree(rpacket_t * packet);

inline double rpacket_get_d_boundfree(rpacket_t * packet);

inline void rpacket_set_d_continuum(rpacket_t * packet, double d_continuum);

inline void rpacket_set_d_electron(rpacket_t * packet, double d_electron);

inline void rpacket_set_d_freefree(rpacket_t * packet, double d_freefree);

inline void rpacket_set_d_boundfree(rpacket_t * packet, double d_boundfree);

inline double rpacket_get_chi_continuum(rpacket_t * packet);

inline double rpacket_get_chi_electron(rpacket_t * packet);

inline double rpacket_get_chi_freefree(rpacket_t * packet);

inline double rpacket_get_chi_boundfree(rpacket_t * packet);

inline void rpacket_set_chi_continuum(rpacket_t * packet, double chi_continuum);

inline void rpacket_set_chi_electron(rpacket_t * packet, double chi_electron);

inline void rpacket_set_chi_freefree(rpacket_t * packet, double chi_freefree);

inline void rpacket_set_chi_boundfree(rpacket_t * packet, double chi_boundfree);

inline double rpacket_get_d_line(rpacket_t * packet);

inline void rpacket_set_d_line(rpacket_t * packet, double d_line);

inline int rpacket_get_next_shell_id(rpacket_t * packet);

inline void rpacket_set_next_shell_id(rpacket_t * packet, int next_shell_id);

inline rpacket_status_t rpacket_get_status(rpacket_t * packet);

inline void rpacket_set_status(rpacket_t * packet, rpacket_status_t status);

inline packet_status_t packet_get_status(rpacket_t * packet);

inline void packet_set_status(rpacket_t * packet, packet_status_t status);

inline int64_t packet_is_r_packet(rpacket_t * packet);

inline int64_t packet_is_k_packet(rpacket_t * packet);

inline int64_t packet_is_i_packet(rpacket_t * packet);

inline tardis_error_t rpacket_init(rpacket_t * packet,
				   storage_model_t * storage, int packet_index,
				   int virtual_packet_flag);

inline void check_array_bounds(int64_t ioned, int64_t nrow, int64_t ncolums);

inline void set_array_int(int64_t irow, int64_t icolums, int64_t nrow,
			  int64_t ncolums, int64_t * array, int64_t val);

inline void set_array_double(int64_t irow, int64_t icolums, int64_t nrow,
			     int64_t ncolums, double *array, double val);

inline int64_t get_array_int(int64_t irow, int64_t icolums, int64_t nrow,
			     int64_t ncolums, int64_t * array);

inline double get_array_double(int64_t irow, int64_t icolums, int64_t nrow,
			       int64_t ncolums, double *array);

inline void calculate_chi_bf(rpacket_t * packet, storage_model_t * storage);

void montecarlo_bound_free_scatter(rpacket_t * packet,
				   storage_model_t * storage, double distance);

void montecarlo_free_free_scatter(rpacket_t * packet, storage_model_t * storage,
				  double distance);

inline montecarlo_event_handler_t montecarlo_continuum_event_handler(rpacket_t *
								     packet,
								     storage_model_t
								     * storage);

inline montecarlo_event_handler_t get_r_event_handler(rpacket_t * packet,
						      storage_model_t * storage,
						      double *distance);

inline montecarlo_event_handler_t get_event_handler(rpacket_t * packet,
						    storage_model_t * storage,
						    double *distance);

inline void rpacket_reset_tau_event(rpacket_t * packet);

inline double rpacket_doppler_factor(rpacket_t * packet,
				     storage_model_t * storage);

inline void rpacket_set_comov_nu(rpacket_t * packet, double comov_nu);

inline double rpacket_get_comov_nu(rpacket_t * packet);

inline double rpacket_get_comov_energy(rpacket_t * packet);

inline void rpacket_set_comov_energy(rpacket_t * packet, double comov_energy);

#endif				// TARDIS_CMONTECARLO_H
