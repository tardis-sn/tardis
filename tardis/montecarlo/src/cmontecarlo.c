#include <inttypes.h>
#ifdef WITHOPENMP
#include <omp.h>
#endif
#include "cmontecarlo.h"

/** Look for a place to insert a value in an inversely sorted float array.
 *
 * @param x an inversely (largest to lowest) sorted float array
 * @param x_insert a value to insert
 * @param imin lower bound
 * @param imax upper bound
 *
 * @return index of the next boundary to the left
 */
static tardis_error_t
reverse_binary_search (const double *x, double x_insert,
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
      int imid = (imin + imax)>>1;
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
	  imid = (imin + imax)>>1;
	}
      if (imax - imin == 2 && x_insert < x[imin + 1])
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

tardis_error_t
line_search (const double *nu, double nu_insert, int64_t number_of_lines,
	     int64_t * result)
{
  tardis_error_t ret_val = TARDIS_ERROR_OK;
  int64_t imin = 0;
  int64_t imax = number_of_lines - 1;
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

static tardis_error_t
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


double
rpacket_doppler_factor (const rpacket_t *packet, const storage_model_t *storage)
{
  return 1.0 -
    rpacket_get_mu (packet) * rpacket_get_r (packet) *
    storage->inverse_time_explosion * INVERSE_C;
}

/* Methods for calculating continuum opacities */
double
bf_cross_section(const storage_model_t * storage, int64_t continuum_id, double comov_nu)
{
  int64_t result;
  tardis_error_t error = binary_search (storage->photo_xsect[continuum_id]->nu, comov_nu, 0,
	       storage->photo_xsect[continuum_id]->no_of_points - 1, &result);
  if (error == TARDIS_ERROR_BOUNDS_ERROR)
    {
      return 0.0;
    }
  else
    {
      double bf_xsect = storage->photo_xsect[continuum_id]->x_sect[result-1]
        + (comov_nu - storage->photo_xsect[continuum_id]->nu[result-1])
        / (storage->photo_xsect[continuum_id]->nu[result] - storage->photo_xsect[continuum_id]->nu[result-1])
        * (storage->photo_xsect[continuum_id]->x_sect[result] - storage->photo_xsect[continuum_id]->x_sect[result-1]);
      return bf_xsect;
    }
}

void calculate_chi_bf(rpacket_t * packet, storage_model_t * storage)
{
  double doppler_factor = rpacket_doppler_factor (packet, storage);
  double comov_nu = rpacket_get_nu (packet) * doppler_factor;

  int64_t no_of_continuum_edges = storage->no_of_edges;
  int64_t current_continuum_id;
  line_search(storage->continuum_list_nu, comov_nu, no_of_continuum_edges, &current_continuum_id);
  rpacket_set_current_continuum_id(packet, current_continuum_id);

  int64_t shell_id = rpacket_get_current_shell_id(packet);
  double T = storage->t_electrons[shell_id];
  double boltzmann_factor = exp(-(H * comov_nu) / (KB*T));

  double bf_helper = 0.0;
  for(int64_t i = current_continuum_id; i < no_of_continuum_edges; i++)
  {
    // get the levelpopulation for the level ijk in the current shell:
    double l_pop = storage->l_pop[shell_id * no_of_continuum_edges + i];
    // get the levelpopulation ratio \frac{n_{0,j+1,k}}{n_{i,j,k}} \frac{n_{i,j,k}}{n_{0,j+1,k}}^{*}:
    double l_pop_r = storage->l_pop_r[shell_id * no_of_continuum_edges + i];
    double bf_x_sect = bf_cross_section(storage, i, comov_nu);
    if (bf_x_sect == 0.0)
      {
        for (int64_t j = i; j < no_of_continuum_edges; j++)
          {
            packet->chi_bf_tmp_partial[j] = bf_helper;
          }
        break;
      }
    bf_helper += l_pop * bf_x_sect * (1.0 - l_pop_r * boltzmann_factor);

    packet->chi_bf_tmp_partial[i] = bf_helper;
  }

  rpacket_set_chi_boundfree(packet, bf_helper * doppler_factor);
}

void calculate_chi_ff(rpacket_t * packet, const storage_model_t * storage)
{
  double doppler_factor = rpacket_doppler_factor (packet, storage);
  double comov_nu = rpacket_get_nu (packet) * doppler_factor;
  int64_t shell_id = rpacket_get_current_shell_id(packet);
  double T = storage->t_electrons[shell_id];
  double boltzmann_factor = exp(-(H * comov_nu) / KB / T);
  double chi_ff_factor = storage->chi_ff_factor[shell_id];

  double chi_ff = chi_ff_factor * (1 - boltzmann_factor) * pow(comov_nu, -3);

  rpacket_set_chi_freefree(packet, chi_ff * doppler_factor);
}

double
compute_distance2boundary (rpacket_t * packet, const storage_model_t * storage)
{
  double r = rpacket_get_r (packet);
  double mu = rpacket_get_mu (packet);
  double r_outer = storage->r_outer[rpacket_get_current_shell_id (packet)];
  double r_inner = storage->r_inner[rpacket_get_current_shell_id (packet)];
  double d_outer =
    sqrt (r_outer * r_outer + ((mu * mu - 1.0) * r * r)) - (r * mu);
  double d_inner;
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

tardis_error_t
compute_distance2line (const rpacket_t * packet, const storage_model_t * storage,
		       double *result)
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
      int64_t cur_zone_id = rpacket_get_current_shell_id (packet);
      double doppler_factor = 1.0 - mu * r * inverse_t_exp * INVERSE_C;
      double comov_nu = nu * doppler_factor;
      if (comov_nu < nu_line)
	{
	  if (rpacket_get_next_line_id (packet) == storage->no_of_lines - 1)
	    {
	      fprintf (stderr, "last_line = %f\n",
		       storage->
		       line_list_nu[rpacket_get_next_line_id (packet) - 1]);
	      fprintf (stderr, "Last line in line list reached!");
	    }
	  else if (rpacket_get_next_line_id (packet) == 0)
	    {
	      fprintf (stderr, "First line in line list!");
	      fprintf (stderr, "next_line = %f\n",
		       storage->
		       line_list_nu[rpacket_get_next_line_id (packet) + 1]);
	    }
	  else
	    {
	      fprintf (stderr, "last_line = %f\n",
		       storage->
		       line_list_nu[rpacket_get_next_line_id (packet) - 1]);
	      fprintf (stderr, "next_line = %f\n",
		       storage->
		       line_list_nu[rpacket_get_next_line_id (packet) + 1]);
	    }
	  fprintf (stderr, "ERROR: Comoving nu less than nu_line!\n");
	  fprintf (stderr, "comov_nu = %f\n", comov_nu);
	  fprintf (stderr, "nu_line = %f\n", nu_line);
	  fprintf (stderr, "(comov_nu - nu_line) / nu_line = %f\n",
		   (comov_nu - nu_line) / nu_line);
	  fprintf (stderr, "r = %f\n", r);
	  fprintf (stderr, "mu = %f\n", mu);
	  fprintf (stderr, "nu = %f\n", nu);
	  fprintf (stderr, "doppler_factor = %f\n", doppler_factor);
	  fprintf (stderr, "cur_zone_id = %" PRIi64 "\n", cur_zone_id);
	  ret_val = TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE;
	}
      else
	{
	  *result = ((comov_nu - nu_line) / nu) * C * t_exp;
	}
    }
  return ret_val;
}

void
compute_distance2continuum(rpacket_t * packet, storage_model_t * storage)
{
  double chi_freefree, chi_electron, chi_continuum, d_continuum;

  if (storage->cont_status == CONTINUUM_ON)
  {
    calculate_chi_bf(packet, storage);
    double chi_boundfree = rpacket_get_chi_boundfree(packet);
    (storage->ff_status == FREE_FREE_ON) ? calculate_chi_ff(packet, storage) : rpacket_set_chi_freefree(packet, 0.0);
    chi_freefree = rpacket_get_chi_freefree(packet);
    chi_electron = storage->electron_densities[packet->current_shell_id] * storage->sigma_thomson *
       rpacket_doppler_factor (packet, storage);
    chi_continuum = chi_boundfree + chi_freefree + chi_electron;
    d_continuum = rpacket_get_tau_event(packet) / chi_continuum;
  }
  else
  {
    chi_electron = storage->electron_densities[rpacket_get_current_shell_id(packet)] * storage->sigma_thomson;
    chi_continuum = chi_electron;
    d_continuum = storage->inverse_electron_densities[rpacket_get_current_shell_id (packet)] *
      storage->inverse_sigma_thomson * rpacket_get_tau_event (packet);
  }

  if (rpacket_get_virtual_packet(packet) > 0)
    {
	  //Set all continuum distances to MISS_DISTANCE in case of an virtual_packet
	  rpacket_set_d_continuum(packet, MISS_DISTANCE);
	  //rpacket_set_chi_boundfree(packet, 0.0);
	  //rpacket_set_chi_electron(packet, chi_electron);
	  //rpacket_set_chi_freefree(packet, 0.0);
      rpacket_set_chi_continuum(packet, chi_continuum);
	}
	else
	{

//        fprintf(stderr, "--------\n");
//        fprintf(stderr, "nu = %e \n", rpacket_get_nu(packet));
//        fprintf(stderr, "chi_electron = %e\n", chi_electron);
//        fprintf(stderr, "chi_boundfree = %e\n", calculate_chi_bf(packet, storage));
//        fprintf(stderr, "chi_line = %e \n", rpacket_get_tau_event(packet) / rpacket_get_d_line(packet));
//        fprintf(stderr, "--------\n");

	  rpacket_set_chi_freefree(packet, chi_freefree);
	  rpacket_set_chi_electron(packet, chi_electron);
	  rpacket_set_chi_continuum(packet, chi_continuum);
	  rpacket_set_d_continuum(packet, d_continuum);
	}
}

int64_t
macro_atom (const rpacket_t * packet, const storage_model_t * storage, rk_state *mt_state)
{
  int emit = 0, i = 0, probability_idx = -1;
  int activate_level =
    storage->line2macro_level_upper[rpacket_get_next_line_id (packet) - 1];
  while (emit != -1)
    {
      double event_random = rk_double (mt_state);
      i = storage->macro_block_references[activate_level] - 1;
      double p = 0.0;
      do
	{

	  probability_idx = ((++i) * storage->no_of_shells +
	  rpacket_get_current_shell_id (packet));
	  p += storage->transition_probabilities[probability_idx];
	}
      while (p <= event_random);
      emit = storage->transition_type[i];
      activate_level = storage->destination_level_id[i];
    }
  return storage->transition_line_id[i];
}

double sample_nu_free_bound(const rpacket_t * packet, const storage_model_t * storage, int64_t continuum_id,
rk_state *mt_state)
{
    double th_frequency = storage->continuum_list_nu[continuum_id];
	int64_t shell_id = rpacket_get_current_shell_id(packet);
	double T = storage->t_electrons[shell_id];
	double zrand = (rk_double(mt_state));
	return th_frequency * (1 - (KB * T / H / th_frequency * log(zrand)));	// Lucy 2003 MC II Eq.26
}

double sample_nu_free_free(const rpacket_t * packet, const storage_model_t * storage, rk_state *mt_state)
{
	int64_t shell_id = rpacket_get_current_shell_id(packet);
	double T = storage->t_electrons[shell_id];
	double zrand = (rk_double(mt_state));
	return -KB * T / H * log(zrand);	// Lucy 2003 MC II Eq.41
}

int64_t
sample_cooling_processes(rpacket_t * packet, rk_state *mt_state, double *individual_cooling_probabilities,
int64_t *cooling_references, int64_t no_of_individual_processes)
{
  int64_t shell_id = rpacket_get_current_shell_id(packet), j = -1;
  double zrand = (rk_double(mt_state)), p = 0.0;

  // TODO: maybe sum before and use bisection
  do
    {
      p += individual_cooling_probabilities[shell_id * no_of_individual_processes + (++j)];
    }
  while(p <= zrand);

  return cooling_references[j];
}

void
k_packet(rpacket_t * packet, const storage_model_t * storage, next_interaction2process * next_process, rk_state *mt_state)
{
  int64_t shell_id = rpacket_get_current_shell_id(packet);
  double zrand = (rk_double(mt_state));
  double cooling_prob = 0.0;
  if (zrand < (cooling_prob += storage->fb_cooling_prob[shell_id]))
    {
      * next_process = BF_EMISSION;
    }
  else if (zrand < (cooling_prob += storage->ff_cooling_prob[shell_id]))
    {
      * next_process = FF_EMISSION;
    }
  else if (zrand < (cooling_prob += storage->coll_exc_cooling_prob[shell_id]))
    {
      * next_process = EXCITATION;
    }
  else
    {
      * next_process = IONIZATION;
    }
  switch(* next_process)
    {
      case FF_EMISSION:
        // No need to select an individual process.
        break;

      case BF_EMISSION:
        {
          int64_t continuum_id = sample_cooling_processes(packet, mt_state, storage->fb_cooling_prob_individual,
            storage->fb_cooling_references, storage->fb_cooling_prob_nd);
          rpacket_set_current_continuum_id(packet, continuum_id);
        }
        break;

      case IONIZATION:
        {
          int64_t activate_level = sample_cooling_processes(packet, mt_state, storage->coll_ion_cooling_prob_individual,
            storage->coll_ion_cooling_references, storage->coll_ion_cooling_prob_nd);
          rpacket_set_macro_atom_activation_level(packet, activate_level);
        }
        break;

      case EXCITATION:
        {
          int64_t activate_level = sample_cooling_processes(packet, mt_state, storage->coll_exc_cooling_prob_individual,
            storage->coll_exc_cooling_references, storage->coll_exc_cooling_prob_nd);
          rpacket_set_macro_atom_activation_level(packet, activate_level);
        }
        break;

      default:
        fprintf(stderr, "This process for k-packet destruction should not exist! (next_process = %d)\n", * next_process);
        exit(1);
    }
}

void
macro_atom_new (rpacket_t * packet, const storage_model_t * storage, next_interaction2process * macro_atom_deactivation_type,
int activation2level_or_cont, rk_state *mt_state)
{
  int emit = activation2level_or_cont;
  int i = 0;
  int64_t activate_level = rpacket_get_macro_atom_activation_level(packet);

  /*
     Do internal jumps until deactivation occurs:
     - radiatively from a macro-atom level (emit = -1)
     - radiatively from a continuum level (emit = -3)
     - collisionally (emit = -2)
  */
  while (emit >= 0)
    {
      double event_random = rk_double (mt_state);
      double p = 0.0;
      if (emit != 2) // Macro-atom is in a normal level.
        {
          int probability_idx = -1;
          i = storage->macro_block_references[activate_level] - 1;
          do
	        {
	          probability_idx = ((++i) * storage->no_of_shells +
	          rpacket_get_current_shell_id (packet));
	          p += storage->transition_probabilities[probability_idx];
	          //fprintf(stderr, " p = %.5e, shell = %d, i=%d \n", p, rpacket_get_current_shell_id (packet), i);
	        }
          while ((p <= event_random));
          emit = storage->transition_type[i];
          activate_level = storage->destination_level_id[i];
          //fprintf(stderr, "MA: activate level:%d, emit:%d \n", activate_level, emit);
        }
      else  // Macro-atom is in a continuum level.
        {
          i = storage->macro_block_references_continuum[activate_level] - 1;
          //fprintf(stderr, "activate_level= %d; j1 = %d\n", activate_level, i);
          //fprintf(stderr, " event_random = %.5e \n", event_random);
          do
	        {
	          p += storage->transition_probabilities_continuum[rpacket_get_current_shell_id (packet) *
				     storage->transition_probabilities_nd_continuum +
				     (++i)];
			  //fprintf(stderr, " p = %.5e, shell = %d, j=%d \n", p, rpacket_get_current_shell_id (packet), i);
	        }
          while ((p <= event_random));
          emit = storage->transition_type_continuum[i];
          activate_level = storage->destination_level_id_continuum[i];
         // fprintf(stderr, "Cont: activate level:%d, emit:%d \n", activate_level, emit);
        }
    }
  switch (emit)
    {
    // radiative deactivation from a level within the macro ion (not a continuum level)
    case -1:
      storage->last_line_interaction_out_id[rpacket_get_id (packet)] = storage->transition_line_id[i];
      * macro_atom_deactivation_type = BB_EMISSION;
      break;

    // radiative deactivation from a continuum level
    case -3:
      // continuum_id of edge corresponding to a continuum transition probability in the macro-atom
      rpacket_set_current_continuum_id(packet, storage->transition_continuum_id[i]);
      * macro_atom_deactivation_type = BF_EMISSION;
      break;

    // collisional deactivation from level or continuum
    case -2:
      * macro_atom_deactivation_type = KPACKET_CREATION;
      break;

    default:
      fprintf(stderr, "This process for macro-atom deactivation should not exist! (emit = %d)\n", emit);
      exit(1);
    }
}


void line_emission(rpacket_t * packet, storage_model_t * storage, rk_state *mt_state)
{
  bool virtual_close_line = false;
  double inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
  int64_t emission_line_id = storage->last_line_interaction_out_id[rpacket_get_id (packet)];
  rpacket_set_nu (packet,
		      storage->line_list_nu[emission_line_id] * inverse_doppler_factor);
  rpacket_set_nu_line (packet, storage->line_list_nu[emission_line_id]);
  rpacket_set_next_line_id (packet, emission_line_id + 1);
  rpacket_reset_tau_event (packet, mt_state);
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
	  montecarlo_one_packet (storage, packet, 1, mt_state);
	  rpacket_set_close_line (packet, old_close_line);
	  virtual_close_line = false;
    }
  test_for_close_line(packet, storage);
}

void bf_emission(rpacket_t * packet, storage_model_t * storage, rk_state *mt_state)
{
  double inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
  int64_t emission_continuum_id = rpacket_get_current_continuum_id(packet);
  double nu_comov = sample_nu_free_bound(packet, storage, emission_continuum_id, mt_state);
  rpacket_set_nu (packet, nu_comov * inverse_doppler_factor);
  rpacket_reset_tau_event (packet, mt_state);
  rpacket_set_recently_crossed_boundary (packet, 0);

  // Have to find current position in line list
  //bool close_line;
  int64_t current_line_id;
  line_search (storage->line_list_nu, nu_comov, storage->no_of_lines, &current_line_id);

  bool last_line = (current_line_id == storage->no_of_lines);
  rpacket_set_next_line_id (packet, current_line_id);
  rpacket_set_last_line (packet, last_line);
  test_for_close_line(packet, storage);
  // Missing: set some interaction ids

  if (rpacket_get_virtual_packet_flag (packet) > 0)
    {
      montecarlo_one_packet (storage, packet, 1, mt_state);
    }
}

void ff_emission(rpacket_t * packet, storage_model_t * storage, rk_state *mt_state)
{
  double inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
  double nu_comov = sample_nu_free_free(packet, storage, mt_state);
  rpacket_set_nu (packet, nu_comov * inverse_doppler_factor);
  rpacket_reset_tau_event (packet, mt_state);
  rpacket_set_recently_crossed_boundary (packet, 0);

  // Have to find current position in line list
  //bool close_line;
  int64_t current_line_id;
  line_search (storage->line_list_nu, nu_comov, storage->no_of_lines, &current_line_id);

  bool last_line = (current_line_id == storage->no_of_lines);
  rpacket_set_next_line_id (packet, current_line_id);
  rpacket_set_last_line (packet, last_line);
  test_for_close_line(packet, storage);
  // Missing: set some interaction ids

  if (rpacket_get_virtual_packet_flag (packet) > 0)
    {
      montecarlo_one_packet (storage, packet, 1, mt_state);
    }
}

void e_packet(rpacket_t * packet, storage_model_t * storage, e_packet_type etype, rk_state *mt_state)
{
  next_interaction2process next_process = etype;

  // Process the e-packet until either bb-, bf- or ff-emission occurs
  while (next_process >= 0)
    {
      switch(next_process)
      {
        case KPACKET_CREATION:
          k_packet(packet, storage, &next_process, mt_state);
          break;

        case EXCITATION:
          macro_atom_new(packet, storage, &next_process, 0, mt_state);
          break;

        case IONIZATION:
          macro_atom_new(packet, storage, &next_process, 2, mt_state);
          break;

        default:
          fprintf (stderr, "Something has gone horribly wrong!\n");
          exit(1);
      }
    }
  // Handle the emission process
  switch (next_process)
   {
     case BB_EMISSION:
       line_emission(packet, storage, mt_state);
       break;

     case BF_EMISSION:
       bf_emission(packet, storage, mt_state);
       break;

     case FF_EMISSION:
       ff_emission(packet, storage, mt_state);
       break;

     default:
       fprintf(stderr, "No emission process was selected.\n");
       exit(1);
   }
}

double
move_packet (rpacket_t * packet, storage_model_t * storage, double distance)
{
  double doppler_factor = rpacket_doppler_factor (packet, storage);
  if (distance > 0.0)
    {
      double r = rpacket_get_r (packet);
      double new_r =
	sqrt (r * r + distance * distance +
	      2.0 * r * distance * rpacket_get_mu (packet));
      rpacket_set_mu (packet,
		      (rpacket_get_mu (packet) * r + distance) / new_r);
      rpacket_set_r (packet, new_r);
      if (rpacket_get_virtual_packet (packet) <= 0)
	{
	  double comov_energy = rpacket_get_energy (packet) * doppler_factor;
	  double comov_nu = rpacket_get_nu (packet) * doppler_factor;
#ifdef WITHOPENMP
#pragma omp atomic
#endif
	  storage->js[rpacket_get_current_shell_id (packet)] +=
	    comov_energy * distance;
#ifdef WITHOPENMP
#pragma omp atomic
#endif
	  storage->nubars[rpacket_get_current_shell_id (packet)] +=
	    comov_energy * distance * comov_nu;
	}
    }
  return doppler_factor;
}

void
increment_j_blue_estimator (const rpacket_t * packet, storage_model_t * storage,
			    double d_line, int64_t j_blue_idx)
{
  double r = rpacket_get_r (packet);
  double r_interaction =
    sqrt (r * r + d_line * d_line +
	  2.0 * r * d_line * rpacket_get_mu (packet));
  double mu_interaction = (rpacket_get_mu (packet) * r + d_line) / r_interaction;
  double doppler_factor = 1.0 - mu_interaction * r_interaction *
    storage->inverse_time_explosion * INVERSE_C;
  double comov_energy = rpacket_get_energy (packet) * doppler_factor;
#ifdef WITHOPENMP
#pragma omp atomic
#endif
  storage->line_lists_j_blues[j_blue_idx] +=
    comov_energy / rpacket_get_nu (packet);
}

int64_t
montecarlo_one_packet (storage_model_t * storage, rpacket_t * packet,
		       int64_t virtual_mode, rk_state *mt_state)
{
  int64_t reabsorbed=-1;
  if (virtual_mode == 0)
    {
      reabsorbed = montecarlo_one_packet_loop (storage, packet, 0, mt_state);
    }
  else
    {
      if ((rpacket_get_nu (packet) > storage->spectrum_virt_start_nu) && (rpacket_get_nu(packet) < storage->spectrum_virt_end_nu))
	{
	  for (int64_t i = 0; i < rpacket_get_virtual_packet_flag (packet); i++)
	    {
              double weight;
              rpacket_t virt_packet = *packet;
              double mu_min;
	      if (rpacket_get_r(&virt_packet) > storage->r_inner[0])
		{
		  mu_min =
		    -1.0 * sqrt (1.0 -
				 (storage->r_inner[0] / rpacket_get_r(&virt_packet)) *
				 (storage->r_inner[0] / rpacket_get_r(&virt_packet)));
		}
	      else
		{
		  mu_min = 0.0;
		}
	      double mu_bin = (1.0 - mu_min) / rpacket_get_virtual_packet_flag (packet);
	      rpacket_set_mu(&virt_packet,mu_min + (i + rk_double (mt_state)) * mu_bin);
	      switch (virtual_mode)
		{
		case -2:
		  weight = 1.0 / rpacket_get_virtual_packet_flag (packet);
		  break;
		case -1:
		  weight =
		    2.0 * rpacket_get_mu(&virt_packet) /
		    rpacket_get_virtual_packet_flag (packet);
		  break;
		case 1:
		  weight =
		    (1.0 -
		     mu_min) / 2.0 / rpacket_get_virtual_packet_flag (packet);
		  break;
		default:
		  fprintf (stderr, "Something has gone horribly wrong!\n");
                  // FIXME MR: we need to somehow signal an error here
                  // I'm adding an exit() here to inform the compiler about the impossible path
                  exit(1);
		}
	      double doppler_factor_ratio =
		rpacket_doppler_factor (packet, storage) /
		rpacket_doppler_factor (&virt_packet, storage);
	      rpacket_set_energy(&virt_packet,
		rpacket_get_energy (packet) * doppler_factor_ratio);
	      rpacket_set_nu(&virt_packet,rpacket_get_nu (packet) * doppler_factor_ratio);
	      reabsorbed = montecarlo_one_packet_loop (storage, &virt_packet, 1, mt_state);
#ifdef WITH_VPACKET_LOGGING
#ifdef WITHOPENMP
#pragma omp critical
		  {
#endif // WITHOPENMP
		    if (storage->virt_packet_count >= storage->virt_array_size)
		      {
			storage->virt_array_size *= 2;
			storage->virt_packet_nus = realloc(storage->virt_packet_nus, sizeof(double) * storage->virt_array_size);
			storage->virt_packet_energies = realloc(storage->virt_packet_energies, sizeof(double) * storage->virt_array_size);
			storage->virt_packet_last_interaction_in_nu = realloc(storage->virt_packet_last_interaction_in_nu, sizeof(double) * storage->virt_array_size);
      storage->virt_packet_last_interaction_type = realloc(storage->virt_packet_last_interaction_type, sizeof(int64_t) * storage->virt_array_size);
      storage->virt_packet_last_line_interaction_in_id = realloc(storage->virt_packet_last_line_interaction_in_id, sizeof(int64_t) * storage->virt_array_size);
      storage->virt_packet_last_line_interaction_out_id = realloc(storage->virt_packet_last_line_interaction_out_id, sizeof(int64_t) * storage->virt_array_size);
		      }
		    storage->virt_packet_nus[storage->virt_packet_count] = rpacket_get_nu(&virt_packet);
		    storage->virt_packet_energies[storage->virt_packet_count] = rpacket_get_energy(&virt_packet) * weight;
        storage->virt_packet_last_interaction_in_nu[storage->virt_packet_count] = storage->last_interaction_in_nu[rpacket_get_id (packet)];
        storage->virt_packet_last_interaction_type[storage->virt_packet_count] = storage->last_interaction_type[rpacket_get_id (packet)];
        storage->virt_packet_last_line_interaction_in_id[storage->virt_packet_count] = storage->last_line_interaction_in_id[rpacket_get_id (packet)];
        storage->virt_packet_last_line_interaction_out_id[storage->virt_packet_count] = storage->last_line_interaction_out_id[rpacket_get_id (packet)];
		    storage->virt_packet_count += 1;
#ifdef WITHOPENMP
		  }
#endif // WITHOPENMP
#endif // WITH_VPACKET_LOGGING
      if ((rpacket_get_nu(&virt_packet) < storage->spectrum_end_nu) &&
          (rpacket_get_nu(&virt_packet) > storage->spectrum_start_nu))
      {
#ifdef WITHOPENMP
#pragma omp critical
		  {
#endif // WITHOPENMP
        int64_t virt_id_nu =
          floor ((rpacket_get_nu(&virt_packet) -
                storage->spectrum_start_nu) /
              storage->spectrum_delta_nu);
        storage->spectrum_virt_nu[virt_id_nu] +=
          rpacket_get_energy(&virt_packet) * weight;
#ifdef WITHOPENMP
		  }
#endif // WITHOPENMP
      }
	    }
	}
      else
	{
	  return 1;
	}
    }
  return reabsorbed;
}

void
move_packet_across_shell_boundary (rpacket_t * packet,
				   storage_model_t * storage, double distance, rk_state *mt_state)
{
  move_packet (packet, storage, distance);
  if (rpacket_get_virtual_packet (packet) > 0)
    {
      double delta_tau_event = rpacket_get_chi_continuum(packet) * distance;
      rpacket_set_tau_event (packet,
			     rpacket_get_tau_event (packet) +
			     delta_tau_event);
    }
  else
    {
      rpacket_reset_tau_event (packet, mt_state);
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
	   (rk_double (mt_state) > storage->inner_boundary_albedo))
    {
      rpacket_set_status (packet, TARDIS_PACKET_STATUS_REABSORBED);
    }
  else
    {
      double doppler_factor = rpacket_doppler_factor (packet, storage);
      double comov_nu = rpacket_get_nu (packet) * doppler_factor;
      double comov_energy = rpacket_get_energy (packet) * doppler_factor;
      rpacket_set_mu (packet, rk_double (mt_state));
      double inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
      rpacket_set_nu (packet, comov_nu * inverse_doppler_factor);
      rpacket_set_energy (packet, comov_energy * inverse_doppler_factor);
      rpacket_set_recently_crossed_boundary (packet, 1);
      if (rpacket_get_virtual_packet_flag (packet) > 0)
	{
	  montecarlo_one_packet (storage, packet, -2, mt_state);
	}
    }
}

void
montecarlo_thomson_scatter (rpacket_t * packet, storage_model_t * storage,
			    double distance, rk_state *mt_state)
{
  double doppler_factor = move_packet (packet, storage, distance);
  double comov_nu = rpacket_get_nu (packet) * doppler_factor;
  double comov_energy = rpacket_get_energy (packet) * doppler_factor;
  rpacket_set_mu (packet, 2.0 * rk_double (mt_state) - 1.0);
  double inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
  rpacket_set_nu (packet, comov_nu * inverse_doppler_factor);
  rpacket_set_energy (packet, comov_energy * inverse_doppler_factor);
  rpacket_reset_tau_event (packet, mt_state);
  rpacket_set_recently_crossed_boundary (packet, 0);
  storage->last_interaction_type[rpacket_get_id (packet)] = 1;
  if (rpacket_get_virtual_packet_flag (packet) > 0)
    {
      montecarlo_one_packet (storage, packet, 1, mt_state);
    }
}

void
montecarlo_bound_free_scatter (rpacket_t * packet, storage_model_t * storage, double distance, rk_state *mt_state)
{
  /* current position in list of continuum edges -> indicates which bound-free processes are possible */
  int64_t current_continuum_id = rpacket_get_current_continuum_id(packet);

  // Determine in which continuum the bf-absorption occurs
  double nu = rpacket_get_nu(packet);
  double chi_bf = rpacket_get_chi_boundfree(packet);
  // get new zrand
  double zrand = rk_double(mt_state);
  double zrand_x_chibf = zrand * chi_bf;

  int64_t ccontinuum = current_continuum_id; /* continuum_id of the continuum in which bf-absorption occurs */

  while ((ccontinuum < storage->no_of_edges - 1) && (packet->chi_bf_tmp_partial[ccontinuum] <= zrand_x_chibf))
  {
    ccontinuum++;
  }
  rpacket_set_current_continuum_id(packet, ccontinuum);
//  Alternative way to choose a continuum for bf-absorption:
//  error =
//  binary_search(storage->chi_bf_tmp_partial, zrand_x_chibf, current_continuum_id,no_of_continuum_edges-1,&ccontinuum);
//  if (error == TARDIS_ERROR_BOUNDS_ERROR) // x_insert < x[imin] -> set index equal to imin
//   {
//      ccontinuum = current_continuum_id;
//   }

  /* Move the packet to the place of absorption, select a direction for re-emission and impose energy conservation
     in the co-moving frame. */
  double old_doppler_factor = move_packet (packet, storage, distance);
  rpacket_set_mu (packet, 2.0 * rk_double (mt_state) - 1.0);
  double inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
  double comov_energy = rpacket_get_energy (packet) * old_doppler_factor;
  rpacket_set_energy (packet, comov_energy * inverse_doppler_factor);
  storage->last_interaction_type[rpacket_get_id (packet)] = 3; // last interaction was a bf-absorption

  // TODO: maybe only do that if we create IONIZATION_ENERGY
  rpacket_set_macro_atom_activation_level(packet,
   storage->cont_edge2macro_continuum[rpacket_get_current_continuum_id(packet)]);

  // Convert the rpacket to thermal or ionization energy
  zrand = (rk_double(mt_state));
  (zrand < storage->continuum_list_nu[ccontinuum] / (nu * rpacket_doppler_factor (packet, storage))) ?
    e_packet(packet, storage, IONIZATION_ENERGY, mt_state): e_packet(packet, storage, THERMAL_ENERGY, mt_state);
}

void
montecarlo_free_free_scatter(rpacket_t * packet, storage_model_t * storage, double distance, rk_state *mt_state)
{
  /* Move the packet to the place of absorption, select a direction for re-emission and impose energy conservation
     in the co-moving frame. */
  double old_doppler_factor = move_packet (packet, storage, distance);
  rpacket_set_mu (packet, 2.0 * rk_double (mt_state) - 1.0);
  double inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
  double comov_energy = rpacket_get_energy (packet) * old_doppler_factor;
  rpacket_set_energy (packet, comov_energy * inverse_doppler_factor);
  storage->last_interaction_type[rpacket_get_id (packet)] = 4; // last interaction was a ff-absorption

  // Create a kpacket
  e_packet(packet, storage, THERMAL_ENERGY, mt_state);
}

void test_for_close_line(rpacket_t * packet, const storage_model_t * storage)
{
  if (!rpacket_get_last_line (packet) &&
      fabs (storage->line_list_nu[rpacket_get_next_line_id (packet)] -
  	    rpacket_get_nu_line (packet)) / rpacket_get_nu_line (packet) <
      1e-7)
    {
      rpacket_set_close_line (packet, true);
    }
}

#ifndef WITH_CONTINUUM
void
montecarlo_line_scatter (rpacket_t * packet, storage_model_t * storage,
			 double distance, rk_state *mt_state)
{
  int64_t line2d_idx = rpacket_get_next_line_id (packet)
  * storage->no_of_shells + rpacket_get_current_shell_id (packet);
  if (rpacket_get_virtual_packet (packet) == 0)
    {
      increment_j_blue_estimator (packet, storage, distance, line2d_idx);
    }
  double tau_line =
    storage->line_lists_tau_sobolevs[line2d_idx];
  double tau_continuum = rpacket_get_chi_continuum(packet) * distance;
  double tau_combined = tau_line + tau_continuum;
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
      double old_doppler_factor = move_packet (packet, storage, distance);
      rpacket_set_mu (packet, 2.0 * rk_double (mt_state) - 1.0);
      double inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
      double comov_energy = rpacket_get_energy (packet) * old_doppler_factor;
      rpacket_set_energy (packet, comov_energy * inverse_doppler_factor);
      storage->last_interaction_in_nu[rpacket_get_id (packet)] =
  rpacket_get_nu (packet);
      storage->last_line_interaction_in_id[rpacket_get_id (packet)] =
	rpacket_get_next_line_id (packet) - 1;
      storage->last_line_interaction_shell_id[rpacket_get_id (packet)] =
	rpacket_get_current_shell_id (packet);
      storage->last_interaction_type[rpacket_get_id (packet)] = 2;
      int64_t emission_line_id = 0;
      if (storage->line_interaction_id == 0)
	{
	  emission_line_id = rpacket_get_next_line_id (packet) - 1;
	}
      else if (storage->line_interaction_id >= 1)
	{
	  emission_line_id = macro_atom (packet, storage, mt_state);
	}
      storage->last_line_interaction_out_id[rpacket_get_id (packet)] =
	emission_line_id;
      rpacket_set_nu (packet,
		      storage->line_list_nu[emission_line_id] *
		      inverse_doppler_factor);
      rpacket_set_nu_line (packet, storage->line_list_nu[emission_line_id]);
      rpacket_set_next_line_id (packet, emission_line_id + 1);
      rpacket_reset_tau_event (packet, mt_state);
      rpacket_set_recently_crossed_boundary (packet, 0);
      if (rpacket_get_virtual_packet_flag (packet) > 0)
	{
	  bool virtual_close_line = false;
	  if (!rpacket_get_last_line (packet) &&
	      fabs (storage->line_list_nu[rpacket_get_next_line_id (packet)] -
		    rpacket_get_nu_line (packet)) <
	      (rpacket_get_nu_line (packet)* 1e-7))
	    {
	      virtual_close_line = true;
	    }
	  // QUESTIONABLE!!!
	  bool old_close_line = rpacket_get_close_line (packet);
	  rpacket_set_close_line (packet, virtual_close_line);
	  montecarlo_one_packet (storage, packet, 1, mt_state);
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
	    rpacket_get_nu_line (packet)) < (rpacket_get_nu_line (packet)*
      1e-7))
    {
      rpacket_set_close_line (packet, true);
    }
}
#else
void
montecarlo_line_scatter (rpacket_t * packet, storage_model_t * storage,
			 double distance, rk_state *mt_state)
{
  int64_t line2d_idx = rpacket_get_next_line_id (packet)
  * storage->no_of_shells + rpacket_get_current_shell_id (packet);
  if (rpacket_get_virtual_packet (packet) == 0)
    {
      increment_j_blue_estimator (packet, storage, distance, line2d_idx);
    }
  double tau_line =
    storage->line_lists_tau_sobolevs[line2d_idx];
  double tau_continuum = rpacket_get_chi_continuum(packet) * distance;
  double tau_combined = tau_line + tau_continuum;
  rpacket_set_next_line_id (packet, rpacket_get_next_line_id (packet) + 1);

  if (rpacket_get_next_line_id (packet) == storage->no_of_lines)
    {
      rpacket_set_last_line (packet, true);
    }
  if (rpacket_get_virtual_packet (packet) > 0)
    {
      rpacket_set_tau_event (packet,
			     rpacket_get_tau_event (packet) + tau_line);
	  test_for_close_line(packet, storage);
    }
  else if (rpacket_get_tau_event (packet) < tau_combined)
    {
      double old_doppler_factor = move_packet (packet, storage, distance);
      rpacket_set_mu (packet, 2.0 * rk_double (mt_state) - 1.0);
      double inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
      double comov_energy = rpacket_get_energy (packet) * old_doppler_factor;
      rpacket_set_energy (packet, comov_energy * inverse_doppler_factor);
      storage->last_interaction_in_nu[rpacket_get_id (packet)] =
  rpacket_get_nu (packet);
      storage->last_line_interaction_in_id[rpacket_get_id (packet)] =
	rpacket_get_next_line_id (packet) - 1;
      storage->last_line_interaction_shell_id[rpacket_get_id (packet)] =
	rpacket_get_current_shell_id (packet);
      storage->last_interaction_type[rpacket_get_id (packet)] = 2;
      int64_t emission_line_id = 0;
      if (storage->line_interaction_id == 0)
	{
	  emission_line_id = rpacket_get_next_line_id (packet) - 1;
	  storage->last_line_interaction_out_id[rpacket_get_id (packet)] =
	    emission_line_id;
	  line_emission(packet, storage, mt_state);
	}
      else if (storage->line_interaction_id >= 1)
	{
	  rpacket_set_macro_atom_activation_level(packet,
	    storage->line2macro_level_upper[rpacket_get_next_line_id (packet) - 1]);
	  e_packet (packet, storage, EXCITATION_ENERGY, mt_state);
	}
	}
  else
    {
      rpacket_set_tau_event (packet,
			     rpacket_get_tau_event (packet) - tau_line);
	  test_for_close_line(packet, storage);
    }
}
#endif // ifndef WITH_CONTINUUM

static void
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
      // FIXME MR: return status of compute_distance2line() is ignored
      rpacket_set_d_line (packet, d_line);
      compute_distance2continuum (packet, storage);
    }
}

static montecarlo_event_handler_t
get_event_handler (rpacket_t * packet, storage_model_t * storage,
		   double *distance, rk_state *mt_state)
{
  montecarlo_compute_distances (packet, storage);
  double d_boundary = rpacket_get_d_boundary (packet);
  double d_continuum = rpacket_get_d_continuum (packet);
  double d_line = rpacket_get_d_line (packet);
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
      handler = montecarlo_continuum_event_handler(packet, storage, mt_state);
    }
  return handler;
}

montecarlo_event_handler_t
montecarlo_continuum_event_handler(rpacket_t * packet, storage_model_t * storage, rk_state *mt_state)
{
  if (storage->cont_status == CONTINUUM_OFF)
    {
      return &montecarlo_thomson_scatter;
    }
  else
    {
  double zrand = (rk_double(mt_state));
  double normaliz_cont_th = rpacket_get_chi_electron(packet)/rpacket_get_chi_continuum(packet);
  double normaliz_cont_bf = rpacket_get_chi_boundfree(packet)/rpacket_get_chi_continuum(packet);

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
}

int64_t
montecarlo_one_packet_loop (storage_model_t * storage, rpacket_t * packet,
			    int64_t virtual_packet, rk_state *mt_state)
{
  rpacket_set_tau_event (packet, 0.0);
  rpacket_set_nu_line (packet, 0.0);
  rpacket_set_virtual_packet (packet, virtual_packet);
  rpacket_set_status (packet, TARDIS_PACKET_STATUS_IN_PROCESS);
  // Initializing tau_event if it's a real packet.
  if (virtual_packet == 0)
    {
      rpacket_reset_tau_event (packet,mt_state);
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
      get_event_handler (packet, storage, &distance, mt_state) (packet, storage,
						      distance, mt_state);
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

void
montecarlo_main_loop(storage_model_t * storage, int64_t virtual_packet_flag, int nthreads, unsigned long seed)
{
  storage->virt_packet_count = 0;
#ifdef WITH_VPACKET_LOGGING
  storage->virt_packet_nus = (double *)malloc(sizeof(double) * storage->no_of_packets);
  storage->virt_packet_energies = (double *)malloc(sizeof(double) * storage->no_of_packets);
  storage->virt_packet_last_interaction_in_nu = (double *)malloc(sizeof(double) * storage->no_of_packets);
  storage->virt_packet_last_interaction_type = (int64_t *)malloc(sizeof(int64_t) * storage->no_of_packets);
  storage->virt_packet_last_line_interaction_in_id = (int64_t *)malloc(sizeof(int64_t) * storage->no_of_packets);
  storage->virt_packet_last_line_interaction_out_id = (int64_t *)malloc(sizeof(int64_t) * storage->no_of_packets);
  storage->virt_array_size = storage->no_of_packets;
#endif // WITH_VPACKET_LOGGING
#ifdef WITHOPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(nthreads);
#pragma omp parallel
  {
    rk_state mt_state;
    rk_seed (seed + omp_get_thread_num(), &mt_state);
    #pragma omp master
      {
        fprintf(stderr, "Running with OpenMP - %d threads\n", omp_get_num_threads());
      }
    #ifdef WITH_CONTINUUM
    double *chi_bf_tmp_partial = calloc(storage->no_of_edges, sizeof(double));
    if (chi_bf_tmp_partial == NULL)
      {
        fprintf(stderr, "Out of memory\n.");
        exit(1);
      }
    #endif // WITH_CONTINUUM

#pragma omp for
#else
  fprintf(stderr, "Running without OpenMP\n");
  rk_state mt_state;
  rk_seed (seed, &mt_state);
  #ifdef WITH_CONTINUUM
    double *chi_bf_tmp_partial = calloc(storage->no_of_edges, sizeof(double));
    if (chi_bf_tmp_partial == NULL)
      {
        fprintf(stderr, "Out of memory\n.");
        exit(1);
      }
  #endif // WITH_CONTINUUM

#endif // WITH_OPENMP
  for (int64_t packet_index = 0; packet_index < storage->no_of_packets; packet_index++)
    {
      int reabsorbed = 0;
      rpacket_t packet;
      rpacket_set_id(&packet, packet_index);
      #ifdef WITH_CONTINUUM
      rpacket_init(&packet, storage, packet_index, virtual_packet_flag, chi_bf_tmp_partial);
      #else
      rpacket_init(&packet, storage, packet_index, virtual_packet_flag);
      #endif // WITH_CONTINUUM
      if (virtual_packet_flag > 0)
	{
	  reabsorbed = montecarlo_one_packet(storage, &packet, -1, &mt_state);
	}
      reabsorbed = montecarlo_one_packet(storage, &packet, 0, &mt_state);
      storage->output_nus[packet_index] = rpacket_get_nu(&packet);
      if (reabsorbed == 1)
	{
	  storage->output_energies[packet_index] = -rpacket_get_energy(&packet);
	}
      else
	{
	  storage->output_energies[packet_index] = rpacket_get_energy(&packet);
	}
    }
  #ifdef WITH_CONTINUUM
  free(chi_bf_tmp_partial);
  #endif // WITH_CONTINUUM
#ifdef WITHOPENMP
  }
#endif
}
