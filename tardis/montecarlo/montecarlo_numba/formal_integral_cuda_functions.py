import numpy as np
from numba import cuda


#Credit for this computation is https://github.com/numba/numba/blob/3fd158f79a12ac5276bc5a72c2404464487c91f0/numba/np/arraymath.py#L3542
@cuda.jit(device=True)
def cuda_searchsorted_value_right(arr, val):
    """
    Find indicies where elements should be inserted 
    to maintain order

    Find the indices into a sorted array a such that, 
    if the corresponding elements in v were inserted 
    before the indices on the right, the order of a 
    would be preserved.

    Parameters
    ----------
    arr : 1-D array_like
    val : 
        Value to be inserted into arr

    Returns
    -------
    int
        location of suitable index
    """
    n = len(arr)
    lo = 0
    hi = n
    while hi > lo:
        mid = (lo + hi) >> 1
        if arr[mid] <= val:
            #mid is too low of an index, go higher
            lo = mid + 1
        else:
            #mid is too high of an index, go down some
            hi = mid
    return lo

@cuda.jit
def cuda_searchsorted(arr, values, output):
    """
    Parameters
    ----------
    arr : 1-D array-like
    values : 1-D array-like
    output : 1-D array-like
        Results of the function being tested given the arr and values
    """
    
    x = cuda.grid(1)
    value = values[x]
    result = cuda_searchsorted_value_right(arr, value)
    cuda.atomic.add(output, x, result)
 

@cuda.jit(device=True)
def cuda_trapz_simpler(arr, dx):
    """
    Parameters
    ----------
    arr1 : (array(float64, 1d, C)
    dx : np.float64
    """
    
    result = arr[0] + arr[-1]
    
    for x in range(1, len(arr) -1):
        result += arr[x] * 2.0
    
    return result * (dx/2.0)

@cuda.jit(debug=True)
def cuda_caller(arr, dx, output):
    """
    arr : (array(float64, 1d, C)
    dx : np.float64
    output : (arr(float64, 1d, C))
        1-d element array
    """
    x = cuda.grid(1)
    output[x] = cuda_trapz_simpler(arr, dx)