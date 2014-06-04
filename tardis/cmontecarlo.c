#include "cmontecarlo.h"

npy_int64 line_search(npy_float64 *nu, npy_float64 nu_insert, npy_int64 number_of_lines)
{
  npy_int64 imin, imax;
  imin = 0;
  imax = number_of_lines - 1;
  if (nu_insert > nu[imin])
    {
      return imin;
    }
  else if (nu_insert < nu[imax])
    {
      return imax + 1;
    }
  else
    {
      return binary_search(nu, nu_insert, imin, imax) + 1;
    }
}

npy_int64 binary_search(npy_float64 *x, npy_float64 x_insert, npy_int64 imin, npy_int64 imax)
{
  if (x_insert > x[imin] || x_insert < x[imax])
    {
      PyErr_SetString(PyExc_ValueError, "Binary Search called but not inside domain. Abort!");
      return -1;
    }
  int imid;
  while (imax - imin > 2)
    {
      imid = (imin + imax) / 2;
      if (x[imid] < x_insert)
	{
	  imax = imid + 1;
	}
      else
	{
	  imin = imid;
	}
    }
  if (imax - imid == 2)
    {
      if (x_insert < x[imin + 1])
	{
	  return imin + 1;
	}
    }
  return imin;
}

inline npy_float64 compute_distance2outer(npy_float64 r, npy_float64 mu, npy_float64 r_outer)
{
  return sqrt(r_outer * r_outer + ((mu * mu - 1.0) * r * r)) - (r * mu);
}

inline npy_float64 compute_distance2inner(npy_float64 r, npy_float64 mu, npy_float64 r_inner)
{
  npy_float64 check = r_inner * r_inner + (r * r * (mu * mu - 1.0));
  if (check < 0.0) 
    {
      return MISS_DISTANCE;
    }
  else
    {
      return mu < 0.0 ? -r * mu - sqrt(check) : MISS_DISTANCE;
    }
}

inline npy_float64 compute_distance2line(npy_float64 r, npy_float64 mu, npy_float64 nu, npy_float64 nu_line, npy_float64 t_exp, npy_float64 inverse_t_exp, npy_float64 last_line, npy_float64 next_line, npy_int64 cur_zone_id)
{
  npy_float64 comov_nu, doppler_factor;
  doppler_factor = 1.0 - mu * r * inverse_t_exp * INVERSE_C;
  comov_nu = nu * doppler_factor;
  if (comov_nu < nu_line)
    {
      fprintf(stderr, "ERROR: Comoving nu less than nu_line!\n");
      fprintf(stderr, "comov_nu = %f\n", comov_nu);
      fprintf(stderr, "nu_line = %f\n", nu_line);
      fprintf(stderr, "(comov_nu - nu_line) / nu_line = %f\n", (comov_nu - nu_line) / nu_line);
      fprintf(stderr, "last_line = %f\n", last_line);
      fprintf(stderr, "next_line = %f\n", next_line);
      fprintf(stderr, "r = %f\n", r);
      fprintf(stderr, "mu = %f\n", mu);
      fprintf(stderr, "nu = %f\n", nu);
      fprintf(stderr, "doppler_factor = %f\n", doppler_factor);
      fprintf(stderr, "cur_zone_id = %d\n", cur_zone_id);
      PyErr_SetString(PyExc_RuntimeError, "comov_nu less than nu_line");
      return 0;
    }
  return ((comov_nu - nu_line) / nu) * C * t_exp;
}

inline npy_float64 compute_distance2electron(npy_float64 r, npy_float64 mu, npy_float64 tau_event, npy_float64 inverse_ne)
{
  return tau_event * inverse_ne;
}

inline npy_int64 macro_atom(npy_int64 activate_level, npy_float64 *p_transition, npy_int64 p_transition_nd, npy_int64 *type_transition, npy_int64 *target_level_id, npy_int64 *target_line_id, npy_int64 *unroll_reference, npy_int64 cur_zone_id)
{
  npy_int64 emit, i = 0;
  npy_float64 p, event_random = 0.0;
  while (1)
    {
      event_random = rk_double(&mt_state);
      i = unroll_reference[activate_level];
      p = 0.0;
      while (1)
	{
	  p = p + p_transition[cur_zone_id * p_transition_nd + i];
	  if (p > event_random)
	    {
	      emit = type_transition[i];
	      activate_level = target_level_id[i];
	      break;
	    }
	  i += 1;
	}
      if (emit == -1)
	{
	  return target_line_id[i];
	}
    }
}

inline npy_float64 move_packet(npy_float64 *r, npy_float64 *mu, npy_float64 nu, npy_float64 energy, npy_float64 distance, npy_float64 *js, npy_float64 *nubars, npy_float64 inverse_t_exp, npy_int64 cur_zone_id, npy_int64 virtual_packet)
{
  npy_float64 new_r, doppler_factor, comov_energy, comov_nu;
  doppler_factor = 1.0 - mu[0] * r[0] * inverse_t_exp * INVERSE_C;
  if (distance <= 0.0)
    {
      return doppler_factor;
    }
  new_r = sqrt(r[0] * r[0] + distance * distance + 2.0 * r[0] * distance * mu[0]);
  mu[0] = (mu[0] * r[0] + distance) / new_r;
  r[0] = new_r;
  if (virtual_packet > 0)
    {
      return doppler_factor;
    }
  comov_energy = energy * doppler_factor;
  comov_nu = nu * doppler_factor;
  js[cur_zone_id] += comov_energy * distance;
  nubars[cur_zone_id] += comov_energy * distance * comov_nu;
  return doppler_factor;
}

inline void increment_j_blue_estimator(npy_int64 *current_line_id, npy_float64 *current_nu, npy_float64 *current_energy, npy_float64 *mu, npy_float64 *r, npy_float64 d_line, npy_int64 j_blue_idx, npy_float64 inverse_time_explosion, npy_float64 *line_lists_j_blues)
{
  npy_float64 comov_energy, comov_nu, r_interaction, mu_interaction, distance, doppler_factor;
  distance = d_line;
  r_interaction = sqrt(r[0] * r[0] + distance * distance + 2 * r[0] * distance * mu[0]);
  mu_interaction = (mu[0] * r[0] + distance) / r_interaction;
  doppler_factor = 1.0 - mu_interaction * r_interaction * inverse_time_explosion * INVERSE_C;
  comov_energy = current_energy[0] * doppler_factor;
  comov_nu = current_nu[0] * doppler_factor;
  line_lists_j_blues[j_blue_idx] += comov_energy / current_nu[0];
}

void init_storage_model(storage_model_t *storage, 
			npy_float64 *packet_nus,
			npy_float64 *packet_mus,
			npy_float64 *packet_energies,
			npy_float64 *r_inner,
			npy_float64 *line_list_nu,
			npy_float64 *output_nus,
			npy_float64 *output_energies,
			npy_float64 time_explosion,
			npy_float64 spectrum_start_nu,
			npy_float64 spectrum_end_nu,
			npy_float64 spectrum_delta_nu,
			npy_float64 *spectrum_virt_nu,
			npy_float64 sigma_thomson,
			npy_float64 *electron_densities,
			npy_float64 *inverse_electron_densities,
			npy_float64 *js,
			npy_float64 *nubars,
			npy_int64 no_of_lines,
			npy_int64 no_of_packets)
{
  storage->packet_nus = packet_nus;
  storage->packet_mus = packet_mus;
  storage->packet_energies = packet_energies;
  storage->r_inner = r_inner;
  storage->output_nus = output_nus;
  storage->output_energies = output_energies;
  storage->time_explosion = time_explosion;
  storage->inverse_time_explosion = 1.0 / time_explosion;
  storage->spectrum_start_nu = spectrum_start_nu;
  storage->spectrum_end_nu = spectrum_end_nu;
  storage->spectrum_delta_nu = spectrum_delta_nu;
  storage->spectrum_virt_nu = spectrum_virt_nu;
  storage->sigma_thomson = sigma_thomson;
  storage->inverse_sigma_thomson = 1.0 / sigma_thomson;
  storage->electron_densities = electron_densities;
  storage->inverse_electron_densities = inverse_electron_densities;
  storage->js = js;
  storage->nubars = nubars;
  storage->no_of_lines = no_of_lines;
  storage->no_of_packets = no_of_packets;
  storage->current_packet_id = -1;
}
