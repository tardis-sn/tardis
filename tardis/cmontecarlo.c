#include "cmontecarlo.h"

rk_state mt_state;

void
initialize_random_kit (unsigned long seed)
{
  rk_seed (seed, &mt_state);
}

inline tardis_error_t
line_search (double *nu, double nu_insert,
	     int64_t number_of_lines, int64_t * result)
{
  tardis_error_t ret_val = TARDIS_ERROR_OK;
  int64_t imin, imax;
  imin = 0;
  imax = number_of_lines - 1;
  if (nu_insert > nu[imin])
    {
      *result = imin;
    }
  else if (nu_insert < nu[imax])
    {
      *result = imax + 1;
    }
  else
    {
      ret_val = reverse_binary_search (nu, nu_insert, imin, imax, result);
      *result = *result + 1;
    }
  return ret_val;
}

inline tardis_error_t
reverse_binary_search (double *x, double x_insert,
		       int64_t imin, int64_t imax, int64_t * result)
{
  /*
     Have in mind that *x points to a reverse sorted array.
     That is large values will have small indices and small ones
     will have large indices.
   */
  tardis_error_t ret_val = TARDIS_ERROR_OK;
  if (x_insert > x[imin] || x_insert < x[imax])
    {
      ret_val = TARDIS_ERROR_BOUNDS_ERROR;
    }
  else
    {
      int imid = (imin + imax) / 2;
      while (imax - imin > 2)
	{
	  if (x[imid] < x_insert)
	    {
	      imax = imid + 1;
	    }
	  else
	    {
	      imin = imid;
	    }
	  imid = (imin + imax) / 2;
	}
      if (imax - imid == 2 && x_insert < x[imin + 1])
	{
	  *result = imin + 1;
	}
      else
	{
	  *result = imin;
	}
    }
  return ret_val;
}

inline tardis_error_t
binary_search (double *x, double x_insert, int64_t imin,
	       int64_t imax, int64_t * result)
{
  /*
     Have in mind that *x points to a sorted array.
     Like [1,2,3,4,5,...]
   */
  int imid;
  tardis_error_t ret_val = TARDIS_ERROR_OK;
  if (x_insert < x[imin] || x_insert > x[imax])
    {
      ret_val = TARDIS_ERROR_BOUNDS_ERROR;
    }
  else
    {
      while (imax >= imin)
	{
	  imid = (imin + imax) / 2;
	  if (x[imid] == x_insert)
	    {
	      *result = imid;
	      break;
	    }
	  else if (x[imid] < x_insert)
	    {
	      imin = imid + 1;
	    }
	  else
	    {
	      imax = imid - 1;
	    }
	}
      if (imax - imid == 2 && x_insert < x[imin + 1])
	{
	  *result = imin;
	}
      else
	{
	  *result = imin;
	}
    }
  return ret_val;
}

inline void
check_array_bounds (int64_t ioned, int64_t nrow, int64_t ncolums)
{
  if (ioned > ((ncolums + 1) * (nrow + 1)))
    {
      fprintf (stderr, "Array Index Out Of Bounds");
      exit (1);
    }
}

inline void
set_array_int (int64_t irow, int64_t icolums, int64_t nrow,
	       int64_t ncolums, int64_t * array, int64_t val)
{
  int64_t ioned = 0;
  ioned = nrow * icolums + irow;
  check_array_bounds (ioned, nrow, ncolums);
  array[ioned] = val;
}

inline void
set_array_double (int64_t irow, int64_t icolums, int64_t nrow,
		  int64_t ncolums, double *array, double val)
{
  int64_t ioned = 0;
  ioned = nrow * icolums + irow;
  check_array_bounds (ioned, nrow, ncolums);
  array[ioned] = val;
}

inline int64_t
get_array_int (int64_t irow, int64_t icolums, int64_t nrow,
	       int64_t ncolums, int64_t * array)
{
  int64_t ioned = 0;
  ioned = nrow * icolums + irow;
  check_array_bounds (ioned, nrow, ncolums);
  return array[ioned];
}

inline double
get_array_double (int64_t irow, int64_t icolums, int64_t nrow,
		  int64_t ncolums, double *array)
{
  int64_t ioned = 0;
  ioned = nrow * icolums + irow;
  check_array_bounds (ioned, nrow, ncolums);
  return array[ioned];
}

inline void
create_kpacket (rpacket_t * packet, storage_model_t * storage)
{
  double comov_energy, comov_nu, doppler_factor;
  doppler_factor = rpacket_doppler_factor (packet, storage);
  comov_energy = rpacket_get_energy (packet) * doppler_factor;
  comov_nu = rpacket_get_nu (packet) * doppler_factor;
  rpacket_set_comov_energy (packet, comov_energy);
  rpacket_set_comov_nu (packet, comov_nu);
  packet->packet_status = TARDIS_K_PACKET_IN_PROCESS;

//ToDo: compute comiving E, mu and, nu and set these values
}

inline void
create_ipacket (rpacket_t * packet, storage_model_t * storage)
{
  double comov_energy, comov_nu, doppler_factor;
  doppler_factor = rpacket_doppler_factor (packet, storage);
  comov_energy = rpacket_get_energy (packet) * doppler_factor;
  comov_nu = rpacket_get_nu (packet) * doppler_factor;
  rpacket_set_comov_energy (packet, comov_energy);
  rpacket_set_comov_nu (packet, comov_nu);
  packet->packet_status = TARDIS_I_PACKET_IN_PROCESS;

//ToDo: compute comiving E, mu and, nu and set these values
}

inline void
create_rpacket (rpacket_t * packet, double comov_nu,
		double comov_mu, storage_model_t * storage)
{
  double current_energy, current_nu, current_mu, inverse_doppler_factor;
  bool last_line;
  int64_t current_line_id;

  tardis_error_t error;
  inverse_doppler_factor =
    1.0 / (1.0 -
	   comov_mu * rpacket_get_r (packet) *
	   storage->inverse_time_explosion * INVERSE_C);
  current_energy = rpacket_get_comov_energy (packet) * inverse_doppler_factor;
  current_nu = comov_nu * inverse_doppler_factor;
  current_mu = comov_mu;

  rpacket_set_energy (packet, current_energy);
  rpacket_set_nu (packet, current_nu);
  rpacket_set_mu (packet, comov_mu);
  rpacket_set_comov_energy (packet, -1.0);	//comoving energy should never be used in an r packet
  rpacket_set_comov_nu (packet, -1.0);	//comoving nu should never be used in an r packet

  error =
    line_search (storage->line_list_nu, comov_nu, storage->no_of_lines,
		 &current_line_id);
  last_line = (current_line_id == storage->no_of_lines);

  rpacket_set_next_line_id (packet, current_line_id);
  rpacket_set_last_line (packet, last_line);
  rpacket_set_recently_crossed_boundary (packet, false);

  packet->packet_status = TARDIS_R_PACKET_IN_PROCESS;

//ToDo: get restframe mu and nu compute comoving
}

inline double
rpacket_doppler_factor (rpacket_t * packet, storage_model_t * storage)
{
  return 1.0 -
    rpacket_get_mu (packet) * rpacket_get_r (packet) *
    storage->inverse_time_explosion * INVERSE_C;
}

inline double
compute_distance2boundary (rpacket_t * packet, storage_model_t * storage)
{
  double r = rpacket_get_r (packet);
  double mu = rpacket_get_mu (packet);
  double r_outer = storage->r_outer[rpacket_get_current_shell_id (packet)];
  double r_inner = storage->r_inner[rpacket_get_current_shell_id (packet)];
  double d_outer =
    sqrt (r_outer * r_outer + ((mu * mu - 1.0) * r * r)) - (r * mu);
  double d_inner;
  double result;
  if (rpacket_get_recently_crossed_boundary (packet) == 1)
    {
      rpacket_set_next_shell_id (packet, 1);
      return d_outer;
    }
  else
    {
      double check = r_inner * r_inner + (r * r * (mu * mu - 1.0));
      if (check < 0.0)
	{
	  rpacket_set_next_shell_id (packet, 1);
	  return d_outer;
	}
      else
	{
	  d_inner = mu < 0.0 ? -r * mu - sqrt (check) : MISS_DISTANCE;
	}
    }
  if (d_inner < d_outer)
    {
      rpacket_set_next_shell_id (packet, -1);
      return d_inner;
    }
  else
    {
      rpacket_set_next_shell_id (packet, 1);
      return d_outer;
    }
}

extern inline tardis_error_t
compute_distance2line (rpacket_t * packet,
		       storage_model_t * storage, double *result)
{
  tardis_error_t ret_val = TARDIS_ERROR_OK;
  if (rpacket_get_last_line (packet))
    {
      *result = MISS_DISTANCE;
    }
  else
    {
      double r = rpacket_get_r (packet);
      double mu = rpacket_get_mu (packet);
      double nu = rpacket_get_nu (packet);
      double nu_line = rpacket_get_nu_line (packet);
      double t_exp = storage->time_explosion;
      double inverse_t_exp = storage->inverse_time_explosion;
      double last_line =
	storage->line_list_nu[rpacket_get_next_line_id (packet) - 1];
      double next_line =
	storage->line_list_nu[rpacket_get_next_line_id (packet) + 1];
      int64_t cur_zone_id = rpacket_get_current_shell_id (packet);
      double comov_nu, doppler_factor;
      doppler_factor = 1.0 - mu * r * inverse_t_exp * INVERSE_C;
      comov_nu = nu * doppler_factor;
      if (comov_nu < nu_line)
	{
	  fprintf (stderr, "ERROR: Comoving nu less than nu_line!\n");
	  fprintf (stderr, "comov_nu = %f\n", comov_nu);
	  fprintf (stderr, "nu_line = %f\n", nu_line);
	  fprintf (stderr, "(comov_nu - nu_line) / nu_line = %f\n",
		   (comov_nu - nu_line) / nu_line);
	  fprintf (stderr, "r = %f\n", r);
	  fprintf (stderr, "mu = %f\n", mu);
	  fprintf (stderr, "nu = %f\n", nu);
	  fprintf (stderr, "doppler_factor = %f\n", doppler_factor);
	  fprintf (stderr, "cur_zone_id = %d\n", cur_zone_id);
	  ret_val = TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE;
	}
      else
	{
	  *result = ((comov_nu - nu_line) / nu) * C * t_exp;
	}
    }
  return ret_val;
}

inline double
compute_distance2electron (rpacket_t * packet, storage_model_t * storage)
{
  if (rpacket_get_virtual_packet (packet) > 0)
    {
      return MISS_DISTANCE;
    }
  double inverse_ne =
    storage->
    inverse_electron_densities[rpacket_get_current_shell_id (packet)] *
    storage->inverse_sigma_thomson;
  return rpacket_get_tau_event (packet) * inverse_ne;
}

inline void
compute_distance2continuum (rpacket_t * packet, storage_model_t * storage)
{
  double chi_boundfree, chi_freefree, chi_electron, chi_continuum,
    d_continuum;
  if (packet->virtual_packet > 0)
    {
      //Set all continuum distances to MISS_DISTANCE in case of an virtual_packet
      rpacket_set_d_continuum (packet, MISS_DISTANCE);
      rpacket_set_chi_boundfree (packet, 0.0);
      rpacket_set_chi_electron (packet, 0.0);
      rpacket_set_chi_freefree (packet, 0.0);
      rpacket_set_chi_continuum (packet, 0.0);
    }
  else
    {
      // Compute the continuum oddities for a real packet
      calculate_chi_bf (packet, storage);

      chi_boundfree = 0.0;	//calculate_chi_bf(packet, storage);
      chi_boundfree = rpacket_get_chi_boundfree (packet);	// For Debug;
      chi_freefree = 0.0;
      chi_electron = storage->electron_densities[packet->current_shell_id] * storage->sigma_thomson;	// For Debugging set * to /
      chi_continuum = chi_boundfree + chi_freefree + chi_electron;
      d_continuum = rpacket_get_tau_event (packet) / chi_continuum;

//        fprintf(stderr, "--------\n");
//        fprintf(stderr, "nu = %e \n", rpacket_get_nu(packet));
//        fprintf(stderr, "chi_electron = %e\n", chi_electron);
//        fprintf(stderr, "chi_boundfree = %e\n", calculate_chi_bf(packet, storage));
//        fprintf(stderr, "chi_line = %e \n", rpacket_get_tau_event(packet) / rpacket_get_d_line(packet));
//        fprintf(stderr, "--------\n");

      rpacket_set_chi_freefree (packet, chi_freefree);
      rpacket_set_chi_electron (packet, chi_electron);
      rpacket_set_d_continuum (packet, d_continuum);
    }
}

inline void
calculate_chi_bf (rpacket_t * packet, storage_model_t * storage)
{
  double bf_helper = 0;
  double chi_bf = 0;
  double nu_th;
  double l_pop_r;
  double l_pop;
  double level_chi;
  double T;
  double nu;
  double delta_nu;
  int64_t last_shell_id;
  int64_t shell_id;
  int64_t i = 0;
  int64_t I;
  int64_t atom;
  int64_t ion;
  int64_t level;

  shell_id = rpacket_get_current_shell_id (packet);
  last_shell_id = rpacket_get_chi_bf_tmp_nu (packet);
  nu = rpacket_get_nu (packet);
  delta_nu = abs (nu - rpacket_get_chi_bf_tmp_nu (packet));

  //check if we can use the stored values
  if (delta_nu > BF_DELTA_NU_RECOMPUTE || shell_id != last_shell_id)
    {
      // set identification values for the tmp data storage
      T = storage->t_electrons[packet->current_shell_id];
      I = storage->chi_bf_index_to_level_nrow;	// This is equal to the number of levels
      rpacket_set_chi_bf_tmp_nu (packet);

      for (i = 0; i <= I; ++i)
	{
	  nu_th = storage->bound_free_th_frequency[i];
	  if (nu_th < nu)
	    {

	      // get the levelpopulation for the level ijk in the current shell
	      l_pop =
		get_array_double (i,
				  packet->current_shell_id,
				  storage->bf_level_population_nrow,
				  storage->bf_level_population_ncolum,
				  storage->bf_level_population);

	      //get the levelpopulation ratio \frac{n_{0,j+1,k}}{n_{i,j,k}} \frac{n_{i,j,k}}{n_{0,j+1,k}}^{*}
	      l_pop_r =
		get_array_double (i,
				  packet->current_shell_id,
				  storage->bf_lpopulation_ratio_nlte_lte_nrow,
				  storage->
				  bf_lpopulation_ratio_nlte_lte_ncolum,
				  storage->bf_lpopulation_ratio_nlte_lte);

	      level_chi =
		l_pop * storage->bf_cross_sections[i] *
		pow ((nu_th / nu),
		     3) * (1 - (l_pop_r * exp (-(H * nu) / KB / T)));
	      bf_helper += level_chi;
	      // DEBUG
	      //            fprintf(stderr, ">>> \n");
	      //            fprintf(stderr, "exp  = %e \n" ,exp(-(H * nu)/KB /T));
	      //            fprintf(stderr, "lpop = %e \n", l_pop);
	      //            fprintf(stderr, "lpop ratio = %e \n", l_pop_r);
	      //            fprintf(stderr, "<<< \n");
	      // END DEBUG
	      //Save the chi bf partial for this packet;  used later for selecting the bf continuum
	    }
	  else
	    {
	      level_chi = 0;
	    }
	  // store the sum of all chis up to the current level. We store the data reversed, this allows us to use the
	  // existing  binary search
	  set_array_double
	    (storage->bf_lpopulation_ratio_nlte_lte_nrow - i, 1,
	     storage->bf_lpopulation_ratio_nlte_lte_nrow, 2,
	     packet->chi_bf_tmp_partial, bf_helper);
	  // store the current chi
	  set_array_double (i,
			    0,
			    storage->bf_lpopulation_ratio_nlte_lte_nrow,
			    2, packet->chi_bf_tmp_partial, level_chi);
	}
      rpacket_set_chi_boundfree (packet, bf_helper);
    }
}

inline int64_t
macro_atom (rpacket_t * packet, storage_model_t * storage)
{
  int emit = 0, i = 0;
  double p, event_random;
  int activate_level =
    storage->line2macro_level_upper[rpacket_get_next_line_id (packet) - 1];
  while (emit != -1)
    {
      event_random = rk_double (&mt_state);
      i = storage->macro_block_references[activate_level] - 1;
      p = 0.0;
      do
	{
	  p +=
	    storage->
	    transition_probabilities[rpacket_get_current_shell_id (packet) *
				     storage->transition_probabilities_nd +
				     (++i)];
	}
      while (p <= event_random);
      emit = storage->transition_type[i];
      activate_level = storage->destination_level_id[i];
    }
  return storage->transition_line_id[i];
}

inline double
move_packet (rpacket_t * packet, storage_model_t * storage, double distance)
{
  double new_r, doppler_factor, comov_energy, comov_nu;
  doppler_factor = rpacket_doppler_factor (packet, storage);
  if (distance > 0.0)
    {
      double r = rpacket_get_r (packet);
      new_r =
	sqrt (r * r + distance * distance +
	      2.0 * r * distance * rpacket_get_mu (packet));
      rpacket_set_mu (packet,
		      (rpacket_get_mu (packet) * r + distance) / new_r);
      rpacket_set_r (packet, new_r);
      if (rpacket_get_virtual_packet (packet) <= 0)
	{
	  comov_energy = rpacket_get_energy (packet) * doppler_factor;
	  comov_nu = rpacket_get_nu (packet) * doppler_factor;
	  storage->js[rpacket_get_current_shell_id (packet)] +=
	    comov_energy * distance;
	  storage->nubars[rpacket_get_current_shell_id (packet)] +=
	    comov_energy * distance * comov_nu;
	}
    }
  return doppler_factor;
}

inline void
increment_j_blue_estimator (rpacket_t * packet,
			    storage_model_t * storage, double d_line,
			    int64_t j_blue_idx)
{
  double comov_energy, r_interaction, mu_interaction, doppler_factor;
  double r = rpacket_get_r (packet);
  r_interaction =
    sqrt (r * r + d_line * d_line +
	  2.0 * r * d_line * rpacket_get_mu (packet));
  mu_interaction = (rpacket_get_mu (packet) * r + d_line) / r_interaction;
  doppler_factor = 1.0 - mu_interaction * r_interaction *
    storage->inverse_time_explosion * INVERSE_C;
  comov_energy = rpacket_get_energy (packet) * doppler_factor;
  storage->line_lists_j_blues[j_blue_idx] +=
    comov_energy / rpacket_get_nu (packet);
}

int64_t
montecarlo_one_packet (storage_model_t * storage, rpacket_t * packet,
		       int64_t virtual_mode)
{
  int64_t i;
  rpacket_t virt_packet;
  double mu_bin;
  double mu_min;
  double doppler_factor_ratio;
  double weight;
  int64_t virt_id_nu;
  int64_t reabsorbed;
  if (virtual_mode == 0)
    {
      reabsorbed = montecarlo_one_packet_loop (storage, packet, 0);
    }
  else
    {
      for (i = 0; i < rpacket_get_virtual_packet_flag (packet); i++)
	{
	  memcpy ((void *) &virt_packet, (void *) packet, sizeof (rpacket_t));
	  if (virt_packet.r > storage->r_inner[0])
	    {
	      mu_min =
		-1.0 * sqrt (1.0 -
			     (storage->r_inner[0] / virt_packet.r) *
			     (storage->r_inner[0] / virt_packet.r));
	    }
	  else
	    {
	      mu_min = 0.0;
	    }
	  mu_bin = (1.0 - mu_min) / rpacket_get_virtual_packet_flag (packet);
	  virt_packet.mu = mu_min + (i + rk_double (&mt_state)) * mu_bin;
	  switch (virtual_mode)
	    {
	    case -2:
	      weight = 1.0 / rpacket_get_virtual_packet_flag (packet);
	      break;
	    case -1:
	      weight =
		2.0 * virt_packet.mu /
		rpacket_get_virtual_packet_flag (packet);
	      break;
	    case 1:
	      weight =
		(1.0 -
		 mu_min) / 2.0 / rpacket_get_virtual_packet_flag (packet);
	      break;
	    default:
	      fprintf (stderr, "Something has gone horribly wrong!\n");
	    }
	  doppler_factor_ratio =
	    rpacket_doppler_factor (packet, storage) /
	    rpacket_doppler_factor (&virt_packet, storage);
	  virt_packet.energy =
	    rpacket_get_energy (packet) * doppler_factor_ratio;
	  virt_packet.nu = rpacket_get_nu (packet) * doppler_factor_ratio;
	  reabsorbed = montecarlo_one_packet_loop (storage, &virt_packet, 1);
	  if ((virt_packet.nu < storage->spectrum_end_nu) &&
	      (virt_packet.nu > storage->spectrum_start_nu))
	    {
	      virt_id_nu =
		floor ((virt_packet.nu -
			storage->spectrum_start_nu) /
		       storage->spectrum_delta_nu);
	      storage->spectrum_virt_nu[virt_id_nu] +=
		virt_packet.energy * weight;
	    }
	}
    }
  return reabsorbed;
}

void
move_packet_across_shell_boundary (rpacket_t * packet,
				   storage_model_t * storage, double distance)
{
  double comov_energy, doppler_factor, comov_nu, inverse_doppler_factor,
    delta_tau_event;
  move_packet (packet, storage, distance);
  if (rpacket_get_virtual_packet (packet) > 0)
    {
      delta_tau_event = rpacket_get_chi_continuum (packet) * distance;
      rpacket_set_tau_event (packet,
			     rpacket_get_tau_event (packet) +
			     delta_tau_event);
    }
  else
    {
      rpacket_reset_tau_event (packet);
    }
  if ((rpacket_get_current_shell_id (packet) < storage->no_of_shells - 1
       && rpacket_get_next_shell_id (packet) == 1)
      || (rpacket_get_current_shell_id (packet) > 0
	  && rpacket_get_next_shell_id (packet) == -1))
    {
      rpacket_set_current_shell_id (packet,
				    rpacket_get_current_shell_id (packet) +
				    rpacket_get_next_shell_id (packet));
      rpacket_set_recently_crossed_boundary (packet,
					     rpacket_get_next_shell_id
					     (packet));
    }
  else if (rpacket_get_next_shell_id (packet) == 1)
    {
      rpacket_set_status (packet, TARDIS_PACKET_STATUS_EMITTED);
    }
  else if ((storage->reflective_inner_boundary == 0) ||
	   (rk_double (&mt_state) > storage->inner_boundary_albedo))
    {
      rpacket_set_status (packet, TARDIS_PACKET_STATUS_REABSORBED);
    }
  else
    {
      doppler_factor = rpacket_doppler_factor (packet, storage);
      comov_nu = rpacket_get_nu (packet) * doppler_factor;
      comov_energy = rpacket_get_energy (packet) * doppler_factor;
      rpacket_set_mu (packet, rk_double (&mt_state));
      inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
      rpacket_set_nu (packet, comov_nu * inverse_doppler_factor);
      rpacket_set_energy (packet, comov_energy * inverse_doppler_factor);
      rpacket_set_recently_crossed_boundary (packet, 1);
      if (rpacket_get_virtual_packet_flag (packet) > 0)
	{
	  montecarlo_one_packet (storage, packet, -2);
	}
    }
}

void
montecarlo_thomson_scatter (rpacket_t * packet, storage_model_t * storage,
			    double distance)
{
  double comov_energy, doppler_factor, comov_nu, inverse_doppler_factor;
  doppler_factor = move_packet (packet, storage, distance);
  comov_nu = rpacket_get_nu (packet) * doppler_factor;
  comov_energy = rpacket_get_energy (packet) * doppler_factor;
  rpacket_set_mu (packet, 2.0 * rk_double (&mt_state) - 1.0);
  inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
  rpacket_set_nu (packet, comov_nu * inverse_doppler_factor);
  rpacket_set_energy (packet, comov_energy * inverse_doppler_factor);
  rpacket_reset_tau_event (packet);
  rpacket_set_recently_crossed_boundary (packet, 0);
  storage->last_interaction_type[storage->current_packet_id] = 1;
  if (rpacket_get_virtual_packet_flag (packet) > 0)
    {
      montecarlo_one_packet (storage, packet, 1);
    }
}

void
montecarlo_line_scatter (rpacket_t * packet, storage_model_t * storage,
			 double distance)
{
  double comov_energy = 0.0;
  int64_t emission_line_id = 0;
  double old_doppler_factor = 0.0;
  double inverse_doppler_factor = 0.0;
  double tau_line = 0.0;
  double tau_continuum = 0.0;
  double tau_combined = 0.0;
  bool virtual_close_line = false;
  int64_t j_blue_idx = -1;
  if (rpacket_get_virtual_packet (packet) == 0)
    {
      j_blue_idx =
	rpacket_get_current_shell_id (packet) *
	storage->line_lists_j_blues_nd + rpacket_get_next_line_id (packet);
      increment_j_blue_estimator (packet, storage, distance, j_blue_idx);
    }
  tau_line =
    storage->line_lists_tau_sobolevs[rpacket_get_current_shell_id
				     (packet) *
				     storage->line_lists_tau_sobolevs_nd
				     + rpacket_get_next_line_id (packet)];
  tau_continuum = rpacket_get_chi_continuum (packet) * distance;
  tau_combined = tau_line + tau_continuum;
  rpacket_set_next_line_id (packet, rpacket_get_next_line_id (packet) + 1);
  if (rpacket_get_next_line_id (packet) == storage->no_of_lines)
    {
      rpacket_set_last_line (packet, true);
    }
  if (rpacket_get_virtual_packet (packet) > 0)
    {
      rpacket_set_tau_event (packet,
			     rpacket_get_tau_event (packet) + tau_line);
    }
  else if (rpacket_get_tau_event (packet) < tau_combined)
    {
      old_doppler_factor = move_packet (packet, storage, distance);
      rpacket_set_mu (packet, 2.0 * rk_double (&mt_state) - 1.0);
      inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
      comov_energy = rpacket_get_energy (packet) * old_doppler_factor;
      rpacket_set_energy (packet, comov_energy * inverse_doppler_factor);
      storage->last_line_interaction_in_id[storage->current_packet_id] =
	rpacket_get_next_line_id (packet) - 1;
      storage->last_line_interaction_shell_id[storage->current_packet_id] =
	rpacket_get_current_shell_id (packet);
      storage->last_interaction_type[storage->current_packet_id] = 2;
      if (storage->line_interaction_id == 0)
	{
	  emission_line_id = rpacket_get_next_line_id (packet) - 1;
	}
      else if (storage->line_interaction_id >= 1)
	{
	  emission_line_id = macro_atom (packet, storage);
	}
      storage->last_line_interaction_out_id[storage->current_packet_id] =
	emission_line_id;
      rpacket_set_nu (packet,
		      storage->line_list_nu[emission_line_id] *
		      inverse_doppler_factor);
      rpacket_set_nu_line (packet, storage->line_list_nu[emission_line_id]);
      rpacket_set_next_line_id (packet, emission_line_id + 1);
      rpacket_reset_tau_event (packet);
      rpacket_set_recently_crossed_boundary (packet, 0);
      if (rpacket_get_virtual_packet_flag (packet) > 0)
	{
	  virtual_close_line = false;
	  if (!rpacket_get_last_line (packet) &&
	      fabs (storage->line_list_nu[rpacket_get_next_line_id (packet)] -
		    rpacket_get_nu_line (packet)) /
	      rpacket_get_nu_line (packet) < 1e-7)
	    {
	      virtual_close_line = true;
	    }
	  // QUESTIONABLE!!!
	  bool old_close_line = rpacket_get_close_line (packet);
	  rpacket_set_close_line (packet, virtual_close_line);
	  montecarlo_one_packet (storage, packet, 1);
	  rpacket_set_close_line (packet, old_close_line);
	  virtual_close_line = false;
	}
    }
  else
    {
      rpacket_set_tau_event (packet,
			     rpacket_get_tau_event (packet) - tau_line);
    }
  if (!rpacket_get_last_line (packet) &&
      fabs (storage->line_list_nu[rpacket_get_next_line_id (packet)] -
	    rpacket_get_nu_line (packet)) / rpacket_get_nu_line (packet) <
      1e-7)
    {
      rpacket_set_close_line (packet, true);
    }
}

void
montecarlo_bound_free_scatter (rpacket_t * packet,
			       storage_model_t * storage, double distance)
{
  tardis_error_t error;

  double nu;
  double nu_bf_edge;
  double nu_edge;
  double chi_bf;
  double chi_sum;
  double zrand;
  double zrand_x_chibf;

  int64_t ccontinuum;
  int64_t i = 0;
  int64_t I;

  //Determine in which continuum the bf-absorption occurs
  nu_edge = rpacket_get_last_bf_edge (packet);
  nu = rpacket_get_nu (packet);
  I = storage->chi_bf_index_to_level_nrow;
  chi_bf = rpacket_get_chi_boundfree (packet);
  // get new zrand
  zrand = (rk_double (&mt_state));
  zrand_x_chibf = zrand * chi_bf;

  error =
    reverse_binary_search (&packet->chi_bf_tmp_partial, zrand_x_chibf, I,
			   2 * I, &ccontinuum);
  // decide whether we go to ionisation energy
  zrand = (rk_double (&mt_state));
  if (zrand < nu / nu_edge)
    {
      // go th the ionisation energy
      // EXTENDE the MACRO ATOM

    }
  else
    {
      //go to the thermal pool
      create_kpacket (packet, storage);
    }

  //decide whether we go to ionisation energy

  rpacket_set_status (packet, TARDIS_PACKET_STATUS_REABSORBED);
}

//////

void
montecarlo_free_free_scatter (rpacket_t * packet, storage_model_t * storage,
			      double distance)
{
  fprintf (stderr,
	   "Ooups, this should not happen! Free free scattering is not implemented yet. Abort!");
  exit (1);
}

extern inline void
montecarlo_compute_distances (rpacket_t * packet, storage_model_t * storage)
{
  // Check if the last line was the same nu as the current line.
  if (rpacket_get_close_line (packet))
    {
      // If so set the distance to the line to 0.0
      rpacket_set_d_line (packet, 0.0);
      // Reset close_line.
      rpacket_set_close_line (packet, false);
    }
  else
    {
      rpacket_set_d_boundary (packet,
			      compute_distance2boundary (packet, storage));
      double d_line;
      compute_distance2line (packet, storage, &d_line);
      rpacket_set_d_line (packet, d_line);
      //ToDo: Remove old rpacket_set_d_electron
      //rpacket_set_d_electron(packet, compute_distance2electron(packet, storage));
      compute_distance2continuum (packet, storage);
    }
}

double
sample_nu_free_free (rpacket_t * packet, storage_model_t * storage)
{
  double T;
  double zrand;
  double shell_id;

  shell_id = rpacket_get_current_shell_id (packet);
  T = storage_get_current_electron_temperature (storage, shell_id);
  zrand = (rk_double (&mt_state));
  return -KB * T / H * log (zrand);	// Lucy 2003 MC II Eq.41
}

double
sample_nu_free_bound (rpacket_t * packet,
		      storage_model_t * storage, int64_t emission_index)
{
  double T;
  double zrand;
  double shell_id;
  double th_frequency;

  //We don't need to know the atom ion or level
  //Look for the corresponding th frequency
  th_frequency = storage->Cr_fb_ijk_th_frequency[emission_index];

  shell_id = rpacket_get_current_shell_id (packet);
  T = storage_get_current_electron_temperature (storage, shell_id);
  zrand = (rk_double (&mt_state));
  return th_frequency * (1 - (KB * T / H / th_frequency * log (zrand)));	//Lucy 2003 MC II Eq.26
}

void
montecarlo_collisional_excitations (rpacket_t
				    * packet,
				    storage_model_t
				    * storage, double distance)
{
  //Goes into an i packet
  exit (1);
}

void
montecarlo_collisional_ionizations (rpacket_t
				    * packet,
				    storage_model_t
				    * storage, double distance)
{
  //Goes into an i packet
  exit (1);
}

void
montecarlo_cooling_free_free (rpacket_t *
			      packet,
			      storage_model_t * storage, double distance)
{
  double comov_nu;
  double current_nu;
  double comov_mu;
  double current_mu;
  double current_r;
  double Cr_ff_jk_max;

  double T;
  double shell_id;

  shell_id = rpacket_get_current_shell_id (packet);
  Cr_ff_jk_max = rpacket_get_Cr_ff_max (packet);
  //Goes into an r packet
  //sample nu (lucy 2003 eq. 41)
  create_rpacket (packet, comov_nu, comov_mu, storage);
  comov_mu = rk_double (&mt_state);
  current_mu = comov_mu;
  current_r = storage_get_r_for_shell (storage, shell_id);
  comov_nu = sample_nu_free_free (packet, storage);
  current_nu =
    comov_nu / (1 -
		(current_mu * current_r *
		 storage->inverse_time_explosion * INVERSE_C));
  rpacket_set_comov_nu (packet, 0);	//Never use the comov nu of an r packet
  rpacket_set_nu (packet, current_nu);
}

void
montecarlo_cooling_free_bound (rpacket_t * packet,
			       storage_model_t * storage, double distance)
{
  tardis_error_t ret_val = TARDIS_ERROR_OK;
  double comov_nu;
  double current_nu;
  double comov_mu;
  double current_mu;
  double current_r;
  double Cr_fb_ijk_max;
  double *Cr_fb_ijk_cumsum_all_current_shell;
  double zrand_normalized;
  double zrand;

  int64_t shell_id;
  int64_t emission_index;
  int64_t emission_level;
  int64_t emission_ion;
  int64_t emission_atom;

  shell_id = rpacket_get_current_shell_id (packet);
  Cr_fb_ijk_max = rpacket_get_Cr_fb_max (packet);

  //Sample bf level
  zrand = (rk_double (&mt_state));
  zrand_normalized = zrand * Cr_fb_ijk_max;

  Cr_fb_ijk_cumsum_all_current_shell =
    storage->Cr_fb_ijk_cumsum_all +
    (storage->Cr_fb_ijk_cumsum_all_nrow * shell_id);

  ret_val =
    binary_search (Cr_fb_ijk_cumsum_all_current_shell, zrand_normalized, 0,
		   storage->Cr_fb_ijk_cumsum_all_nrow, &emission_index);

  //Goes into an r packet
  //sample nu (lucy 2003 eq. 41)
  create_rpacket (packet, comov_nu, comov_mu, storage);

  comov_mu = rk_double (&mt_state);
  current_mu = comov_mu;
  current_r = storage_get_r_for_shell (storage, shell_id);
  comov_nu = sample_nu_free_bound (packet, storage, emission_index);
  current_nu =
    comov_nu / (1 -
		(current_mu * current_r *
		 storage->inverse_time_explosion * INVERSE_C));
  rpacket_set_comov_nu (packet, 0);	//Never use the comov nu of an r packet
  rpacket_set_nu (packet, current_nu);

} inline montecarlo_event_handler_t

get_k_event_handler (rpacket_t * packet,
		     storage_model_t * storage, double *distance)
{
  double Cr_fb_ijk_max;
  double Cr_ff_jk_max;
  double Cr_bb_ijk_max;
  double Cr_ion_ijk_max;
  double Cr_all_max;
  double zrand, zrand_normalized;

  montecarlo_event_handler_t handler;

  // The k packet performs always an cooling process.
  // We have to select the process via sampling the coolingrates
  Cr_fb_ijk_max =
    get_array_double (storage->Cr_fb_ijk_cumsum_all_nrow - 1,
		      packet->current_shell_id,
		      storage->Cr_fb_ijk_cumsum_all_nrow,
		      storage->Cr_fb_ijk_cumsum_all_ncolum,
		      storage->Cr_fb_ijk_cumsum_all);

  Cr_ff_jk_max =
    get_array_double (storage->Cr_ff_jk_cumsum_all_nrow - 1,
		      packet->current_shell_id,
		      storage->Cr_ff_jk_cumsum_all_nrow,
		      storage->Cr_ff_jk_cumsum_all_ncolum,
		      storage->Cr_ff_jk_cumsum_all);

  Cr_bb_ijk_max =
    get_array_double (storage->Cr_bb_ijk_cumsum_all_nrow - 1,
		      packet->current_shell_id,
		      storage->Cr_bb_ijk_cumsum_all_nrow,
		      storage->Cr_bb_ijk_cumsum_all_ncolum,
		      storage->Cr_bb_ijk_cumsum_all);

  Cr_ion_ijk_max =
    get_array_double (storage->Cr_ion_ijk_cumsum_all_nrow - 1,
		      packet->current_shell_id,
		      storage->Cr_ion_ijk_cumsum_all_nrow,
		      storage->Cr_ion_ijk_cumsum_all_ncolum,
		      storage->Cr_ion_ijk_cumsum_all);

  rpacket_set_Cr_fb_max (packet, Cr_fb_ijk_max);
  rpacket_set_Cr_ff_max (packet, Cr_ff_jk_max);
  rpacket_set_Cr_bb_max (packet, Cr_bb_ijk_max);
  rpacket_set_Cr_ion_max (packet, Cr_ion_ijk_max);

  Cr_all_max = Cr_fb_ijk_max + Cr_ff_jk_max + Cr_bb_ijk_max + Cr_ion_ijk_max;

  //sample the cooling process
  zrand = (rk_double (&mt_state));
  zrand_normalized = zrand * Cr_all_max;
  if (zrand_normalized < Cr_bb_ijk_max)
    {
      //Do bound bound (col exi)
    }
  else if (zrand_normalized < (Cr_bb_ijk_max + Cr_ion_ijk_max))
    {
      //Do col ion
    }
  else if (zrand_normalized < (Cr_bb_ijk_max + Cr_ion_ijk_max + Cr_ff_jk_max))
    {
      //Do free-free
      handler = &montecarlo_cooling_free_free;
    }
  else
    {
      //Do free-bound
      handler = &montecarlo_cooling_free_bound;

    }
  return handler;
}

inline montecarlo_event_handler_t
montecarlo_continuum_event_handler (rpacket_t *
				    packet, storage_model_t * storage)
{
  double zrand, normaliz_cont_th, normaliz_cont_bf, normaliz_cont_ff;
  zrand = (rk_double (&mt_state));
  normaliz_cont_th =
    rpacket_get_chi_electron (packet) / rpacket_get_chi_continuum (packet);
  normaliz_cont_bf =
    rpacket_get_chi_boundfree (packet) / rpacket_get_chi_continuum (packet);
  normaliz_cont_ff =
    rpacket_get_chi_freefree (packet) / rpacket_get_chi_continuum (packet);

  if (zrand < normaliz_cont_th)
    {
      //Return the electron scatter event function
      return &montecarlo_thomson_scatter;
    }
  else if (zrand < (normaliz_cont_th + normaliz_cont_bf))
    {
      //Return the bound-free scatter event function
      return &montecarlo_bound_free_scatter;
    }
  else
    {
      //Return the free-free scatter event function
      return &montecarlo_free_free_scatter;
    }

}

inline montecarlo_event_handler_t
get_r_event_handler (rpacket_t * packet,
		     storage_model_t * storage, double *distance)
{
  double d_boundary, d_continuum, d_line;
  montecarlo_compute_distances (packet, storage);
  d_boundary = rpacket_get_d_boundary (packet);
  //d_continuum = rpacket_get_d_continuum(packet);
  d_continuum = 1e99;		//Do not use the new part of the code ToDo: Remove this later
  d_line = rpacket_get_d_line (packet);

  montecarlo_event_handler_t handler;
  if (d_line <= d_boundary && d_line <= d_continuum)
    {
      *distance = d_line;
      handler = &montecarlo_line_scatter;
    }
  else if (d_boundary <= d_continuum)
    {
      *distance = d_boundary;
      handler = &move_packet_across_shell_boundary;
    }
  else
    {
      *distance = d_continuum;
      handler = montecarlo_continuum_event_handler (packet, storage);
    }
  return handler;
}

inline montecarlo_event_handler_t
get_event_handler (rpacket_t * packet,
		   storage_model_t * storage, double *distance)
{
  montecarlo_event_handler_t handler;
  if (packet_is_r_packet (packet))
    {
      handler = get_r_event_handler (packet, storage, distance);
    }
  else if (packet_is_k_packet (packet))
    {
      handler = get_k_event_handler (packet, storage, distance);
    }
  else if (packet_is_i_packet (packet))
    {
      //ToDo:add the i event handler
    }
  return handler;
//ToDo: add the new event handler for the different packet types
}

int64_t
montecarlo_one_packet_loop (storage_model_t * storage, rpacket_t * packet,
			    int64_t virtual_packet)
{
  rpacket_set_tau_event (packet, 0.0);
  rpacket_set_nu_line (packet, 0.0);
  rpacket_set_virtual_packet (packet, virtual_packet);
  rpacket_set_status (packet, TARDIS_PACKET_STATUS_IN_PROCESS);
  // Initializing tau_event if it's a real packet.
  if (virtual_packet == 0)
    {
      rpacket_reset_tau_event (packet);
    }
  // For a virtual packet tau_event is the sum of all the tau's that the packet passes.
  while (rpacket_get_status (packet) == TARDIS_PACKET_STATUS_IN_PROCESS)
    {
      // Check if we are at the end of line list.
      if (!rpacket_get_last_line (packet))
	{
	  rpacket_set_nu_line (packet,
			       storage->
			       line_list_nu[rpacket_get_next_line_id
					    (packet)]);
	}
      double distance;
      get_event_handler (packet, storage, &distance) (packet, storage,
						      distance);
      if (virtual_packet > 0 && rpacket_get_tau_event (packet) > 10.0)
	{
	  rpacket_set_tau_event (packet, 100.0);
	  rpacket_set_status (packet, TARDIS_PACKET_STATUS_EMITTED);
	}
    }
  if (virtual_packet > 0)
    {
      rpacket_set_energy (packet,
			  rpacket_get_energy (packet) * exp (-1.0 *
							     rpacket_get_tau_event
							     (packet)));
    }
  return rpacket_get_status (packet) ==
    TARDIS_PACKET_STATUS_REABSORBED ? 1 : 0;
}


tardis_error_t
rpacket_init (rpacket_t * packet, storage_model_t * storage, int packet_index,
	      int virtual_packet_flag)
{
  double nu_line;
  double current_r;
  double comov_mu;
  double comov_nu;
  double comov_current_nu;
  double comov_energy;
  double current_nu, current_energy, current_mu;
  int64_t current_line_id;
  int current_shell_id;
  bool last_line;
  bool close_line;
  int recently_crossed_boundary;

  //assign temporary opacity storage pointer from storage to packet
  packet->chi_bf_tmp_partial = storage->chi_bf_tmp_partial;
  packet->chi_bf_tmp_partial_last_shell_id = -1;
  packet->chi_bf_tmp_partial_last_nu = 0;
  //int nr,nc;
  //nr = storage->bf_level_population_nrow;
  //packet.chi_bf_partial = (double *) malloc(nt *  sizeof(double));

  /* Get packet information from storage */
  tardis_error_t ret_val = TARDIS_ERROR_OK;
  comov_nu = storage->packet_nus[packet_index];
  comov_energy = storage->packet_energies[packet_index];
  comov_mu = storage->packet_mus[packet_index];
  comov_current_nu = comov_nu;
  current_shell_id = 0;
  current_r = storage->r_inner[0];
  current_mu = comov_mu;

  current_nu =
    comov_nu / (1 -
		(current_mu * current_r *
		 storage->inverse_time_explosion * INVERSE_C));

  current_energy =
    comov_energy / (1 -
		    (current_mu * current_r *
		     storage->inverse_time_explosion * INVERSE_C));

  if ((ret_val =
       line_search (storage->line_list_nu, comov_current_nu,
		    storage->no_of_lines,
		    &current_line_id)) != TARDIS_ERROR_OK)
    {
      return ret_val;
    }
  /*Set the values of the new r packet corresponding to the packet_index */
  //rpacket_set_nu(packet, current_nu);
  //rpacket_set_mu(packet, current_mu);
  //rpacket_set_energy(packet, current_energy);
  rpacket_set_comov_energy (packet, comov_energy);
  rpacket_set_energy (packet, current_energy);
  rpacket_set_r (packet, current_r);
  rpacket_set_current_shell_id (packet, current_shell_id);
  rpacket_set_close_line (packet, false);
  rpacket_set_recently_crossed_boundary (packet, true);
  rpacket_set_virtual_packet_flag (packet, virtual_packet_flag);
  /*Call create new */
  create_rpacket (packet, comov_nu, comov_mu, storage);
  /*
     last_line = (current_line_id == storage->no_of_lines);o
     recently_crossed_boundary = true;o
     rpacket_set_nu(packet, current_nu);o
     rpacket_set_mu(packet, current_mu);o
     rpacket_set_energy(packet, current_energy);o
     rpacket_set_r(packet, current_r);o
     rpacket_set_current_shell_id(packet, current_shell_id);o
     rpacket_set_next_line_id(packet, current_line_id);o
     rpacket_set_last_line(packet, last_line);o
     rpacket_set_close_line(packet, false);o
     rpacket_set_recently_crossed_boundary(packet, recently_crossed_boundary);
     rpacket_set_virtual_packet_flag(packet, virtual_packet_flag);
   */
  return ret_val;
}

/*
  Getter and setter methods.
*/

extern inline double
rpacket_get_nu (rpacket_t * packet)
{
  return packet->nu;
}

inline void
rpacket_set_nu (rpacket_t * packet, double nu)
{
  packet->nu = nu;
} inline double

rpacket_get_mu (rpacket_t * packet)
{

  return packet->mu;
}

inline void
rpacket_set_mu (rpacket_t * packet, double mu)
{
  packet->mu = mu;
} extern inline double

rpacket_get_energy (rpacket_t * packet)
{

  return packet->energy;
}

inline void
rpacket_set_energy (rpacket_t * packet, double energy)
{
  packet->energy = energy;
} inline double

rpacket_get_r (rpacket_t * packet)
{

  return packet->r;
}

inline void
rpacket_set_r (rpacket_t * packet, double r)
{
  packet->r = r;
} inline double

rpacket_get_tau_event (rpacket_t * packet)
{

  return packet->tau_event;
}

inline void
rpacket_set_tau_event (rpacket_t * packet, double tau_event)
{
  packet->tau_event = tau_event;
} inline double

rpacket_get_nu_line (rpacket_t * packet)
{

  return packet->nu_line;
}

inline void
rpacket_set_nu_line (rpacket_t * packet, double nu_line)
{
  packet->nu_line = nu_line;
} inline unsigned int

rpacket_get_current_shell_id (rpacket_t * packet)
{

  return packet->current_shell_id;
}

inline double
storage_get_current_electron_temperature (storage_model_t *
					  storage, int64_t current_shell_id)
{
  return storage->t_electrons[current_shell_id];
}

inline double
storage_get_r_for_shell (storage_model_t * storage, int64_t current_shell_id)
{
  return storage->r_inner[current_shell_id];
}

inline void
rpacket_set_current_shell_id (rpacket_t * packet,
			      unsigned int current_shell_id)
{
  packet->current_shell_id = current_shell_id;
} inline unsigned int

rpacket_get_next_line_id (rpacket_t * packet)
{

  return packet->next_line_id;
}

inline void
rpacket_set_next_line_id (rpacket_t * packet, unsigned int next_line_id)
{
  packet->next_line_id = next_line_id;
} inline bool

rpacket_get_last_line (rpacket_t * packet)
{

  return packet->last_line;
}

inline void
rpacket_set_last_line (rpacket_t * packet, bool last_line)
{
  packet->last_line = last_line;
}

inline bool
rpacket_get_close_line (rpacket_t * packet)
{
  return packet->close_line;
}

inline void
rpacket_set_close_line (rpacket_t * packet, bool close_line)
{
  packet->close_line = close_line;
}

inline int
rpacket_get_recently_crossed_boundary (rpacket_t * packet)
{
  return packet->recently_crossed_boundary;
}

inline void
rpacket_set_recently_crossed_boundary (rpacket_t * packet,
				       int recently_crossed_boundary)
{
  packet->recently_crossed_boundary = recently_crossed_boundary;
} inline int

rpacket_get_virtual_packet_flag (rpacket_t * packet)
{

  return packet->virtual_packet_flag;
}

inline void
rpacket_set_virtual_packet_flag (rpacket_t * packet, int virtual_packet_flag)
{
  packet->virtual_packet_flag = virtual_packet_flag;
} inline int

rpacket_get_virtual_packet (rpacket_t * packet)
{

  return packet->virtual_packet;
}

inline void
rpacket_set_virtual_packet (rpacket_t * packet, int virtual_packet)
{
  packet->virtual_packet = virtual_packet;
} inline double

rpacket_get_d_boundary (rpacket_t * packet)
{

  return packet->d_boundary;
}

inline void
rpacket_set_d_boundary (rpacket_t * packet, double d_boundary)
{
  packet->d_boundary = d_boundary;
} inline double

rpacket_get_d_line (rpacket_t * packet)
{

  return packet->d_line;
}

inline void
rpacket_set_d_line (rpacket_t * packet, double d_line)
{
  packet->d_line = d_line;
} inline int

rpacket_get_next_shell_id (rpacket_t * packet)
{

  return packet->next_shell_id;
}

inline void
rpacket_set_next_shell_id (rpacket_t * packet, int next_shell_id)
{
  packet->next_shell_id = next_shell_id;
} inline double

rpacket_get_d_continuum (rpacket_t * packet)
{

  return packet->d_cont;
}

inline double
rpacket_get_d_electron (rpacket_t * packet)
{
  return packet->d_th;
}

inline double
rpacket_get_d_freefree (rpacket_t * packet)
{
  return packet->d_ff;
}

inline double
rpacket_get_d_boundfree (rpacket_t * packet)
{
  return packet->d_bf;
}

inline void
rpacket_set_d_continuum (rpacket_t * packet, double d_continuum)
{
  packet->d_cont = d_continuum;
} inline void

rpacket_set_d_electron (rpacket_t * packet, double d_electron)
{

  packet->d_th = d_electron;
} inline void

rpacket_set_d_freefree (rpacket_t * packet, double d_freefree)
{

  packet->d_ff = d_freefree;
} inline void

rpacket_set_d_boundfree (rpacket_t * packet, double d_boundfree)
{

  packet->d_bf = d_boundfree;
} inline double

rpacket_get_chi_continuum (rpacket_t * packet)
{

  return packet->chi_cont;
}

inline double
rpacket_get_chi_electron (rpacket_t * packet)
{
  return packet->chi_th;
}

inline double
rpacket_get_chi_freefree (rpacket_t * packet)
{
  return packet->chi_ff;
}

inline double
rpacket_get_chi_boundfree (rpacket_t * packet)
{
  return packet->chi_bf;
}

inline void
rpacket_set_chi_continuum (rpacket_t * packet, double chi_continuum)
{
  packet->chi_cont = chi_continuum;
} inline void

rpacket_set_chi_electron (rpacket_t * packet, double chi_electron)
{

  packet->chi_th = chi_electron;
} inline void

rpacket_set_chi_freefree (rpacket_t * packet, double chi_freefree)
{

  packet->chi_ff = chi_freefree;
} inline void

rpacket_set_chi_boundfree (rpacket_t * packet, double chi_boundfree)
{
  packet->chi_bf = chi_boundfree;
} inline rpacket_status_t

rpacket_get_status (rpacket_t * packet)
{

  return packet->status;
}

inline void
rpacket_set_status (rpacket_t * packet, rpacket_status_t status)
{
  packet->status = status;
}

inline packet_status_t
packet_get_status (rpacket_t * packet)
{
  return packet->packet_status;
}

inline void
packet_set_status (rpacket_t * packet, packet_status_t status)
{
  packet->status = status & (15);	// get the first 4 bits
  packet->packet_status = status;
  /* if (status & 1<<3)// is in process
     {
     packet->status = TARDIS_PACKET_STATUS_IN_PROCESS;
     }
     else if ( status & 1<<2) // is emitted
     {
     packet->status = TARDIS_PACKET_STATUS_EMITTED;
     }
     else // is reabsorbed
     {
     packet->status = TARDIS_PACKET_STATUS_REABSORBED;
     }

     packet->packet_status = status;
   */
}

/* Other accessor methods. */

inline int64_t
packet_is_r_packet (rpacket_t * packet)
{
  return packet_get_status (packet) & 1 << 6;
}

inline int64_t
packet_is_k_packet (rpacket_t * packet)
{
  //return packet_get_status(packet) & 1 << 5;
  return 0;			// Do not use the new part of the code ToDO: Remove this later
}

inline int64_t
packet_is_i_packet (rpacket_t * packet)
{
  //return packet_get_status(packet) & 1 << 4;
  return 0;			// Do not use the new part of the code ToDO: Remove this later
}

/* Other accessor methods. */

inline void
rpacket_reset_tau_event (rpacket_t * packet)
{
  rpacket_set_tau_event (packet, -log (rk_double (&mt_state)));
}

inline void
rpacket_set_comov_nu (rpacket_t * packet, double comov_nu)
{
  packet->comov_nu = comov_nu;
} inline double

rpacket_get_comov_nu (rpacket_t * packet)
{

  return packet->comov_nu;
}

inline double
rpacket_get_comov_energy (rpacket_t * packet)
{
  return packet->comov_energy;
}

inline void
rpacket_set_comov_energy (rpacket_t * packet, double comov_energy)
{
  packet->comov_energy = comov_energy;
}

inline void
rpacket_set_last_bf_edge (rpacket_t * packet, double nu_edge)
{

  packet->last_bf_edge = nu_edge;
}

inline double
rpacket_get_last_bf_edge (rpacket_t * packet)
{

  return packet->last_bf_edge;
}

inline int64_t
rpacket_get_chi_bf_tmp_shell_id (rpacket_t * packet)
{
/*
Get the shell_id for which we have computed the chi_bf_tmp_partial
*/
  return packet->chi_bf_tmp_partial_last_shell_id;
}

inline void
rpacket_set_chi_bf_tmp_shell_id (rpacket_t * packet, int64_t shell_id)
{
/*
Set the shell_id for which we have computed the chi_bf_tmp_partial
*/
  packet->chi_bf_tmp_partial_last_shell_id = shell_id;
}

inline double
rpacket_get_chi_bf_tmp_nu (rpacket_t * packet)
{
/*
Get the nu for which we have computed the chi_bf_tmp_partial
*/
  return packet->chi_bf_tmp_partial_last_nu;
}

inline void
rpacket_set_chi_bf_tmp_last_nu (rpacket_t * packet, double nu)
{
/*
Set the nu for which we have computed the chi_bf_tmp_partial
*/
  packet->chi_bf_tmp_partial_last_nu = nu;
}

inline void
rpacket_set_Cr_fb_max (rpacket_t * packet, double Cr_fb_max)
{
/*
Set the sum of all bound free cooling rates
*/
  packet->Cr_fb_max = Cr_fb_max;
}

inline double
rpacket_get_Cr_fb_max (rpacket_t * packet)
{
/*
Get the sum of all bound free cooling rates
*/
  return packet->Cr_fb_max;
}

inline void
rpacket_set_Cr_ff_max (rpacket_t * packet, double Cr_ff_max)
{
/*
Set the sum of all free free cooling rates
*/
  packet->Cr_ff_max = Cr_ff_max;
}

inline double
rpacket_get_Cr_ff_max (rpacket_t * packet)
{
/*
Get the sum of all free free cooling rates
*/
  return packet->Cr_ff_max;
}

inline void
rpacket_set_Cr_bb_max (rpacket_t * packet, double Cr_bb_max)
{
/*
Set the sum of all bound bound cooling rates (collisional excitation)
*/
  packet->Cr_bb_max = Cr_bb_max;
}

inline double
rpacket_get_Cr_bb_max (rpacket_t * packet)
{
/*
Get the sum of all bound bound cooling rates (collisional excitation)
*/
  return packet->Cr_bb_max;
}

inline void
rpacket_set_Cr_ion_max (rpacket_t * packet, double Cr_ion_max)
{
/*
Set the sum of all collisional ionization rates
*/
  packet->Cr_ion_max = Cr_ion_max;
}

inline double
rpacket_get_Cr_ion_max (rpacket_t * packet)
{
/*
Get the sum of all collisional ionization rates
*/
  return packet->Cr_ion_max;
}
