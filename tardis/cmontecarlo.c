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
  doppler_factor 1.0 - mu_interaction * r_interaction * inverse_time_explosion * INVERSE_C;
  comov_energy = current_energy[0] * doppler_factor;
  comov_nu = current_nu[0] * doppler_factor;
  line_lists_j_blues[j_blue_idx] += comov_energy / current_nu[0];
}
