#ifdef WITHOPENMP
#include <omp.h>
#endif
#include "cmontecarlo.h"

rk_state mt_state;

void
initialize_random_kit (unsigned long seed)
{
  rk_seed (seed, &mt_state);
}

INLINE tardis_error_t
line_search (double *nu, double nu_insert, int64_t number_of_lines,
	     int64_t * result)
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

INLINE double
rpacket_doppler_factor (rpacket_t * packet, storage_model_t * storage)
{
  return 1.0 -
    rpacket_get_mu (packet) * rpacket_get_r (packet) *
    storage->inverse_time_explosion * INVERSE_C;
}

/* Methods for calculating continuum opacities */

INLINE double
bf_cross_section(storage_model_t * storage, int64_t continuum_id, double comov_nu)
{
  /* Temporary hardcoded values */
  double chi_bf_partial = 0.5 * 0.25e-15;
  double cont_chi_bf[] = {chi_bf_partial, 0.0, 2.0 * chi_bf_partial, 0.3 * chi_bf_partial, 2.0 * chi_bf_partial};
  /* End of temporary hardcoded values */
  double sigma_bf = 0.6 * chi_bf_partial;
  //double sigma_bf = cont_chi_bf[continuum_id]; //storage->bf_cross_sections[continuum_id]
  return sigma_bf * pow((storage->continuum_list_nu[continuum_id] / comov_nu), 3);
}

INLINE
void calculate_chi_bf(rpacket_t * packet, storage_model_t * storage)
{
  double bf_helper = 0;
  double comov_nu, doppler_factor;
  double T;
  double boltzmann_factor;
  int64_t shell_id;
  int64_t current_continuum_id;
  int64_t i;
  int64_t no_of_continuum_edges = storage->no_of_edges;

  doppler_factor = rpacket_doppler_factor (packet, storage);
  comov_nu = rpacket_get_nu (packet) * doppler_factor;

  line_search(storage->continuum_list_nu, comov_nu, no_of_continuum_edges, &current_continuum_id);
  rpacket_set_current_continuum_id(packet, current_continuum_id);

  shell_id = rpacket_get_current_shell_id(packet);
  T = storage->t_electrons[shell_id];
  boltzmann_factor = exp(-(H * comov_nu) / KB / T);

  for(i = current_continuum_id; i < no_of_continuum_edges; i++)
  {
    // get the levelpopulation for the level ijk in the current shell:
    double l_pop = storage->l_pop[shell_id * no_of_continuum_edges + i];
    // get the levelpopulation ratio \frac{n_{0,j+1,k}}{n_{i,j,k}} \frac{n_{i,j,k}}{n_{0,j+1,k}}^{*}:
    double l_pop_r = storage->l_pop_r[shell_id * no_of_continuum_edges + i];
    bf_helper += l_pop * bf_cross_section(storage, i, comov_nu) * (1 - l_pop_r * boltzmann_factor);

    storage->chi_bf_tmp_partial[i] = bf_helper;
  }

  rpacket_set_chi_boundfree(packet, bf_helper * doppler_factor);
}

INLINE double
compute_distance2boundary (rpacket_t * packet, storage_model_t * storage)
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

INLINE tardis_error_t
compute_distance2line (rpacket_t * packet, storage_model_t * storage,
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
      double comov_nu, doppler_factor;
      doppler_factor = 1.0 - mu * r * inverse_t_exp * INVERSE_C;
      comov_nu = nu * doppler_factor;
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
	  fprintf (stderr, "cur_zone_id = %lld\n", cur_zone_id);
	  ret_val = TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE;
	}
      else
	{
	  *result = ((comov_nu - nu_line) / nu) * C * t_exp;
	}
    }
  return ret_val;
}

INLINE void
compute_distance2continuum(rpacket_t * packet, storage_model_t * storage)
{
  double chi_boundfree, chi_freefree, chi_electron, chi_continuum, d_continuum;

  if (storage->cont_status == CONTINUUM_ON)
  {
    calculate_chi_bf(packet, storage);
    chi_boundfree = rpacket_get_chi_boundfree(packet);
    rpacket_set_chi_freefree(packet, 0.0);
    chi_freefree = rpacket_get_chi_freefree(packet);
    chi_electron = storage->electron_densities[packet->current_shell_id] * storage->sigma_thomson *
       rpacket_doppler_factor (packet, storage);
    chi_continuum = chi_boundfree + chi_freefree + chi_electron;
    d_continuum = rpacket_get_tau_event(packet) / chi_continuum;
  }
  else
  {
    chi_electron = storage->electron_densities[packet->current_shell_id] * storage->sigma_thomson;
    chi_continuum = chi_electron;
    d_continuum = storage->inverse_electron_densities[rpacket_get_current_shell_id (packet)] *
      storage->inverse_sigma_thomson * rpacket_get_tau_event (packet);
  }

  if (packet->virtual_packet > 0)
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

INLINE int64_t
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



double sample_nu_free_bound(rpacket_t * packet, storage_model_t * storage, int64_t continuum_id)
{
	double T;
	double zrand;
	int64_t shell_id;
	double th_frequency;

	th_frequency = storage->continuum_list_nu[continuum_id];

	shell_id = rpacket_get_current_shell_id(packet);
	T = storage->t_electrons[shell_id];
	zrand = (rk_double(&mt_state));
	return th_frequency * (1 - (KB * T / H / th_frequency * log(zrand)));	// Lucy 2003 MC II Eq.26
}

#if 0
double sample_nu_free_free(rpacket_t * packet, storage_model_t * storage)
{
	double T;
	double zrand;
	int64_t shell_id;

	shell_id = rpacket_get_current_shell_id(packet);
	T = storage->t_electrons[shell_id];
	zrand = (rk_double(&mt_state));
	return -KB * T / H * log(zrand);	// Lucy 2003 MC II Eq.41
}
#endif

INLINE void
macro_atom_new (rpacket_t * packet, storage_model_t * storage, next_interaction2process * macro_atom_deactivation_type,
int activation2level_or_cont)
{
  int level_or_cont = activation2level_or_cont;
  int emit = 0, i = 0, j = 0, activate_level;
  int64_t emission_line_id = 0;
  int64_t emission_continuum_id;
  double p, event_random;

  switch (activation2level_or_cont)
    {
    // Macro-atom is activated to a normal level.
    case 0:
      activate_level =
      storage->line2macro_level_upper[rpacket_get_next_line_id (packet) - 1];
      break;

     // Macro-atom is activated to a continuum level.
    case 1:
      activate_level =
      storage->cont_edge2macro_continuum[rpacket_get_current_continuum_id(packet)]; // ? -1 as in : rpacket_get_next_line_id (packet) - 1
      fprintf(stderr, "cc_id = %d ", rpacket_get_current_continuum_id(packet));
      fprintf(stderr, "activate_level = %d", activate_level);
      break;
    }

  /*
     Do internal jumps until deactivation occurs:
     - radiatively from a macro-atom level (emit = -1)
     - radiatively from a continuum level (emit = -3)
     - collisionally (emit = -2)
  */
  while (emit >= 0)
    {
      event_random = rk_double (&mt_state);
      p = 0.0;
      if (level_or_cont == 0) // Macro-atom is in a normal level.
        {
          i = storage->macro_block_references[activate_level] - 1;
          do
	        {
	          p += storage->transition_probabilities[rpacket_get_current_shell_id (packet) *
				     storage->transition_probabilities_nd +
				     (++i)];
	        }
          while ((p <= event_random));
          emit = storage->transition_type[i];
          activate_level = storage->destination_level_id[i];
          if (emit == 2) // internal jump to higher ionization state
            {
              level_or_cont = 1; // set macro-atom to operate in continuum
              fprintf(stderr, "Internal jumps to higher ionization states are not implemented yet.\n");
            }
        }
      else  // Macro-atom is in a continuum level.
        {
          j = storage->macro_block_references_continuum[activate_level] - 1; // - 1 because of i++
          fprintf(stderr, " j1 = %d ", j);
          do
	        {
	          p += storage->transition_probabilities_continuum[rpacket_get_current_shell_id (packet) *
				     storage->transition_probabilities_nd_continuum +
				     (++j)];
	        }
          while ((p <= event_random));
          fprintf(stderr, " j2 = %d", j);
          emit = storage->transition_type_continuum[j];
          activate_level = storage->destination_level_id_continuum[j];
          level_or_cont = 0; // set macro-atom to normal level
          // for debug
          if (emit >=0) {fprintf(stderr, "-%d-> A ", activate_level);}
        }
    }
  switch (emit)
    {
    // radiative deactivation from a level within the macro ion (not a continuum level)
    case -1:
      emission_line_id  = storage->transition_line_id[i];
      storage->last_line_interaction_out_id[rpacket_get_id (packet)] = emission_line_id;
      * macro_atom_deactivation_type = BB_EMISSION;
      break;

    // radiative deactivation from a continuum level
    case -3:
      // continuum_id of edge corresponding to a continuum transition probability in the macro-atom
      emission_continuum_id = storage->transition_continuum_id[i];
      rpacket_set_current_continuum_id(packet, emission_continuum_id);
      * macro_atom_deactivation_type = BF_EMISSION;
      break;

    // collisional deactivation from level or continuum
    case -2:
      * macro_atom_deactivation_type = KPACKET_CREATION;
      fprintf(stderr, "Collisional macro-atom deactivations are not implemented yet.\n");
      break;

    default:
      fprintf(stderr, "This process for macro-atom deactivation should not exist!\n");
    }
}


INLINE void line_emission(rpacket_t * packet, storage_model_t * storage)
{
  bool virtual_close_line = false;
  double inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
  int64_t emission_line_id = storage->last_line_interaction_out_id[rpacket_get_id (packet)];
  rpacket_set_nu (packet,
		      storage->line_list_nu[emission_line_id] * inverse_doppler_factor);
  rpacket_set_nu_line (packet, storage->line_list_nu[emission_line_id]);
  rpacket_set_next_line_id (packet, emission_line_id + 1);
  rpacket_reset_tau_event (packet);
  rpacket_set_recently_crossed_boundary (packet, 0);

  // for debug
  fprintf(stderr, "-bb %d-> r\n", emission_line_id);

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
  test_for_close_line(packet, storage);
}

INLINE void bf_emission(rpacket_t * packet, storage_model_t * storage)
{
  double inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
  int64_t emission_continuum_id = rpacket_get_current_continuum_id(packet);
  double nu_comov = sample_nu_free_bound(packet, storage, emission_continuum_id);
  rpacket_set_nu (packet, nu_comov * inverse_doppler_factor);
  rpacket_reset_tau_event (packet);
  rpacket_set_recently_crossed_boundary (packet, 0);

  // Have to find current position in line list
  bool last_line;
  bool close_line;
  int64_t current_line_id;
  line_search (storage->line_list_nu, nu_comov,
		    storage->no_of_lines,
		    &current_line_id);
  last_line = (current_line_id == storage->no_of_lines);
  rpacket_set_next_line_id (packet, current_line_id);
  rpacket_set_last_line (packet, last_line);
  rpacket_set_close_line (packet, false); // ? is this the right thing to do
  // Missing: set some interaction ids

  if (rpacket_get_virtual_packet_flag (packet) > 0)
    {
      montecarlo_one_packet (storage, packet, 1);
    }
}

void e_packet(rpacket_t * packet, storage_model_t * storage, e_packet_type etype)
{
  next_interaction2process next_process;
  switch(etype)
  {
    case EXCITATION_ENERGY:
      // Activate macro-atom to a normal level (not continuum)
      fprintf(stderr, "r --> A");
      macro_atom_new(packet, storage, &next_process, 0);
      break;

    case IONIZATION_ENERGY:
      // Activate macro-atom to a continuum level
      fprintf(stderr, "r --> A* ");
      macro_atom_new(packet, storage, &next_process, 1);
      break;

    case THERMAL_ENERGY:
      //create_kpacket(packet, storage, &next_process);
      //break;
      fprintf(stderr, "r --> k --> reabsorbed\n");
      rpacket_set_status (packet, TARDIS_PACKET_STATUS_REABSORBED);
      return;
  }

  // Process the e-packet until either bb-, bf- or ff-emission occurs
  while (next_process >= 0)
    {
      switch(next_process)
      {
        case KPACKET_CREATION:
          //create_kpacket(packet, storage, &next_process);
          fprintf(stderr, " That should not happen. We cannot create kpackets.\n");
          return;

        case COLL_EXCITATION:
          macro_atom_new(packet, storage, &next_process, 0);
          break;

        case COLL_IONIZATION:
          macro_atom_new(packet, storage, &next_process, 1);
          break;
      }
    }
  // Handle the emission process
  switch (next_process)
   {
     case BB_EMISSION:
       line_emission(packet, storage);
       break;

     case BF_EMISSION:
       fprintf(stderr, "-bf-> r\n");
       bf_emission(packet, storage);
       break;

     case FF_EMISSION:
       fprintf(stderr, " Free-free emissions are not implemented yet.\n");
       break;

     default:
       fprintf(stderr, "No emission process was selected.\n");
   }
}

INLINE double
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

INLINE void
increment_j_blue_estimator (rpacket_t * packet, storage_model_t * storage,
			    double d_line, int64_t j_blue_idx)
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
#ifdef WITHOPENMP
#pragma omp atomic
#endif
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
      if ((rpacket_get_nu (packet) > storage->spectrum_virt_start_nu) && (rpacket_get_nu(packet) < storage->spectrum_virt_end_nu))
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
      else
	{
	  return 1;
	}
    }
  return reabsorbed;
}

void
move_packet_across_shell_boundary (rpacket_t * packet,
				   storage_model_t * storage, double distance)
{
  double comov_energy, doppler_factor, comov_nu, inverse_doppler_factor;
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
  storage->last_interaction_type[rpacket_get_id (packet)] = 1;
  if (rpacket_get_virtual_packet_flag (packet) > 0)
    {
      montecarlo_one_packet (storage, packet, 1);
    }
}

void
montecarlo_bound_free_scatter (rpacket_t * packet, storage_model_t * storage, double distance)
{
  /* current position in list of continuum edges -> indicates which bound-free processes are possible */
  int64_t current_continuum_id = rpacket_get_current_continuum_id(packet);
  int64_t ccontinuum; /* continuum_id of the continuum in which bf-absorption occurs */

  double zrand, zrand_x_chibf, chi_bf, nu;
  // Determine in which continuum the bf-absorption occurs
  nu = rpacket_get_nu(packet); // frequency from before moving the packet

  chi_bf = rpacket_get_chi_boundfree(packet);
  // get new zrand
  zrand = (rk_double(&mt_state));
  zrand_x_chibf = zrand * chi_bf;

  ccontinuum = current_continuum_id;
  fprintf (stderr, "start selecting continuum at %d\n", ccontinuum);
  while (storage->chi_bf_tmp_partial[ccontinuum] <= zrand_x_chibf)
  {
    ccontinuum++;
  }
  fprintf (stderr, "end selecting continuum at %d\n", ccontinuum);
  rpacket_set_current_continuum_id(packet, ccontinuum);
//  Alternative way to choose a continuum for bf-absorption:
//  error =
//  binary_search(storage->chi_bf_tmp_partial, zrand_x_chibf, current_continuum_id,no_of_continuum_edges-1,&ccontinuum);
//  if (error == TARDIS_ERROR_BOUNDS_ERROR) // x_insert < x[imin] -> set index equal to imin
//   {
//      ccontinuum = current_continuum_id;
//   }

  // Move the packet to the place of absorption and impose energy conservation in the co-moving frame.
  double old_doppler_factor;
  double inverse_doppler_factor;
  double comov_energy;
  old_doppler_factor = move_packet (packet, storage, distance);
  rpacket_set_mu (packet, 2.0 * rk_double (&mt_state) - 1.0);
  inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
  comov_energy = rpacket_get_energy (packet) * old_doppler_factor;
  rpacket_set_energy (packet, comov_energy * inverse_doppler_factor);
  storage->last_interaction_type[rpacket_get_id (packet)] = 3; // last interaction was a bf-absorption

  // Convert the rpacket to thermal or ionization energy
  zrand = (rk_double(&mt_state));
  (zrand < storage->continuum_list_nu[ccontinuum] / (nu * rpacket_doppler_factor (packet, storage))) ?
    e_packet(packet, storage, IONIZATION_ENERGY): e_packet(packet, storage, THERMAL_ENERGY);
}

void
montecarlo_free_free_scatter(rpacket_t * packet, storage_model_t * storage, double distance)
{
  rpacket_set_status (packet, TARDIS_PACKET_STATUS_REABSORBED);
}

INLINE void test_for_close_line(rpacket_t * packet, storage_model_t * storage)
{
  if (!rpacket_get_last_line (packet) &&
      fabs (storage->line_list_nu[rpacket_get_next_line_id (packet)] -
  	    rpacket_get_nu_line (packet)) / rpacket_get_nu_line (packet) <
      1e-7)
    {
      rpacket_set_close_line (packet, true);
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
  int64_t j_blue_idx = -1;
  if (rpacket_get_virtual_packet (packet) == 0)
    {
      j_blue_idx =
	rpacket_get_current_shell_id (packet) *
	storage->line_lists_j_blues_nd + rpacket_get_next_line_id (packet);
      increment_j_blue_estimator (packet, storage, distance, j_blue_idx);
    }
  tau_line =
    storage->line_lists_tau_sobolevs[rpacket_get_current_shell_id (packet) *
				     storage->line_lists_tau_sobolevs_nd +
				     rpacket_get_next_line_id (packet)];
  tau_continuum = rpacket_get_chi_continuum(packet) * distance;
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
	  test_for_close_line(packet, storage);
    }
  else if (rpacket_get_tau_event (packet) < tau_combined)
    {
      old_doppler_factor = move_packet (packet, storage, distance);
      rpacket_set_mu (packet, 2.0 * rk_double (&mt_state) - 1.0);
      inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
      comov_energy = rpacket_get_energy (packet) * old_doppler_factor;
      rpacket_set_energy (packet, comov_energy * inverse_doppler_factor);
      storage->last_line_interaction_in_id[rpacket_get_id (packet)] =
	rpacket_get_next_line_id (packet) - 1;
      storage->last_line_interaction_shell_id[rpacket_get_id (packet)] =
	rpacket_get_current_shell_id (packet);
      storage->last_interaction_type[rpacket_get_id (packet)] = 2;
      if (storage->line_interaction_id == 0)
	{
	  emission_line_id = rpacket_get_next_line_id (packet) - 1;
	  storage->last_line_interaction_out_id[rpacket_get_id (packet)] =
	    emission_line_id;
	  line_emission(packet, storage);
	}
      else if (storage->line_interaction_id >= 1)
	{
	  e_packet (packet, storage, EXCITATION_ENERGY);
	}
	}
  else
    {
      rpacket_set_tau_event (packet,
			     rpacket_get_tau_event (packet) - tau_line);
	  test_for_close_line(packet, storage);
    }
}

INLINE void
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
      compute_distance2continuum (packet, storage);
    }
}

INLINE montecarlo_event_handler_t
get_event_handler (rpacket_t * packet, storage_model_t * storage,
		   double *distance)
{
  double d_boundary, d_continuum, d_line;
  montecarlo_compute_distances (packet, storage);
  d_boundary = rpacket_get_d_boundary (packet);
  d_continuum = rpacket_get_d_continuum (packet);
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
      handler = montecarlo_continuum_event_handler(packet, storage);
    }
  return handler;
}

INLINE montecarlo_event_handler_t
montecarlo_continuum_event_handler(rpacket_t * packet, storage_model_t * storage)
{
  if (storage->cont_status == CONTINUUM_OFF)
    {
      return &montecarlo_thomson_scatter;
    }
  else
    {
  double zrand, normaliz_cont_th, normaliz_cont_bf, normaliz_cont_ff;
  zrand = (rk_double(&mt_state));
  normaliz_cont_th = rpacket_get_chi_electron(packet)/rpacket_get_chi_continuum(packet);
  normaliz_cont_bf = rpacket_get_chi_boundfree(packet)/rpacket_get_chi_continuum(packet);
  normaliz_cont_ff = rpacket_get_chi_freefree(packet)/rpacket_get_chi_continuum(packet);

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

void
montecarlo_main_loop(storage_model_t * storage, int64_t virtual_packet_flag, int nthreads, unsigned long seed)
{
  int64_t packet_index;
#ifdef WITHOPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(nthreads);
#pragma omp parallel
  {
    initialize_random_kit(seed + omp_get_thread_num());
#pragma omp for
#else
  initialize_random_kit(seed);
#endif
  for (packet_index = 0; packet_index < storage->no_of_packets; packet_index++)
    {
      int reabsorbed = 0;
      rpacket_t packet;
      rpacket_set_id(&packet, packet_index);
      rpacket_init(&packet, storage, packet_index, virtual_packet_flag);
      if (virtual_packet_flag > 0)
	{
	  reabsorbed = montecarlo_one_packet(storage, &packet, -1);
	}
      reabsorbed = montecarlo_one_packet(storage, &packet, 0);
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
#ifdef WITHOPENMP
  }
#endif
}
