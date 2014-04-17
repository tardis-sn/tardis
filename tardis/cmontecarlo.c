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
