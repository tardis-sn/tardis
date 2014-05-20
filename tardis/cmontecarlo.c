#include "cmontecarlo.h"

npy_int64 binary_search(npy_float64 *x, npy_float64 x_insert, npy_int64 imin, npy_int64 imax)
{
  if (x_insert > x[imin] || x_insert < x[imax])
    {
      PyErr_SetString(PyExc_ValueError, "Binary Search called but not inside domain. Abort!");
      // Not possible to raise an exception, I think.
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
