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

npy_float64 compute_distance2outer(npy_float64 r, npy_float64 mu, npy_float64 r_outer)
{
  return sqrt(r_outer * r_outer + ((mu * mu - 1.0) * r * r)) - (r * mu);
}

npy_float64 compute_distance2inner(npy_float64 r, npy_float64 mu, npy_float64 r_inner)
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

npy_float64 compute_distance2line(npy_float64 r, npy_float64 mu, npy_float64 nu, npy_float64 nu_line, npy_float64 t_exp, npy_float64 inverse_t_exp, npy_float64 last_line, npy_float64 next_line, npy_int64 cur_zone_id)
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

npy_float64 compute_distance2electron(npy_float64 r, npy_float64 mu, npy_float64 tau_event, npy_float64 inverse_ne)
{
  return tau_event * inverse_ne;
}
