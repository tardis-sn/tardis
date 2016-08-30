"""
Unit tests for methods in `tardis/montecarlo/src/integrator.c`.
* `ctypes` library is used to wrap C methods and expose them to python.


Probable Reasons for Failing Tests:
-----------------------------------

1. ???

2. Return type of any method changed:
  - Modify the `restype` parameter in the test method here.
  - For example:
        ```
        cmontecarlo_methods.rpacket_doppler_factor.restype = c_double
        ```

3. Underlying logic modified:
  - Check whether the changes made in C method are logically correct.
  - If the changes made were correct and necessary, update the corresponding
    test case.


General Test Design Procedure:
------------------------------

Please follow this design procedure while adding a new test:

1. Parametrization as per desire of code coverage.
  - C tests have different flows controlled by conditional statements.
    Parameters checked in conditions can be provided in different testcases.
  - Keep consistency with variable names as (in order):
    - `packet_params`
    - `model_params`
    - `expected_params` (`expected` if only one value to be asserted.)
  - Suggested variable names can be compromised if readability of the test
    increases.

2. Test Method body:
  - Keep name as `test_` + `(name of C method)`.
  - Refer to method `test_rpacket_doppler_factor` below for description.
"""

import os
import pytest
from ctypes import CDLL, byref, c_uint, c_int, c_int64, c_double, c_ulong, POINTER
from numpy.testing import assert_equal, assert_almost_equal

from tardis import __path__ as path

cmontecarlo_filepath = os.path.join(path[0], 'montecarlo', 'montecarlo.so')
cmontecarlo_methods = CDLL(cmontecarlo_filepath)

@pytest.mark.parametrize(
    ['cr_idx', 'no_of_cr_shells', 'expected'],
    [( 0  ,  1   ,  -1),
     ( 1  ,  1   ,   1),
     ( 2  ,  3   ,  -1),
     ( 0  ,  3   ,  -1),
     ( 5  ,  3   ,   1)]
)
def test_get_cr_sign(cr_idx, no_of_cr_shells,expected):
    cmontecarlo_methods.get_cr_sign.restype = c_int
    result = cmontecarlo_methods.get_cr_sign(c_int64(cr_idx), c_int64(no_of_cr_shells))
    assert_equal(result,expected)


@pytest.mark.parametrize(
    ['no_of_cr_shells', 'p' , 'R_ph', 'expected'],
    [( 3 , 1797120000000000.0 , 1235520000000000.0, 0),
     ( 3 ,  898560000000000.0 , 1235520000000000.0, 3),
     ( 3 , 1235520000000000.0 , 1235520000000000.0, 0)]
)
def test_get_cr_start(no_of_cr_shells, p, R_ph, expected):
    cmontecarlo_methods.get_cr_start.restype = c_int64
    result = cmontecarlo_methods.get_cr_start(c_int64(no_of_cr_shells), c_double(p), c_double(R_ph))
    assert_equal(result,expected)

@pytest.mark.parametrize(
    ['cr_idx', 'no_of_cr_shells', 'expected'],
    [( 0  ,  1  ,  0),
     ( 1  ,  1  ,  0), # Special case, used when indexing Rs
     ( 2  ,  4  ,  2),
     ( 3  ,  4  ,  3),
     ( 4  ,  4  ,  2),
     ( 6  ,  4  ,  0),
     ( 4  ,  4  ,  2),
     ( 3  ,  3  ,  1),
     ( 2  ,  3  ,  2)]
)
def test_get_sh_idx(cr_idx, no_of_cr_shells,expected):
    cmontecarlo_methods.get_sh_idx.restype = c_int64
    result = cmontecarlo_methods.get_sh_idx(c_int64(cr_idx), c_int64(no_of_cr_shells))
    assert_equal(result,expected)

@pytest.mark.parametrize(
    ['len', 'p' , 'Rs', 'expected'],
    [( 5 , 0.2 , (c_double * 5)(*[0.61, 0.5, 0.3, 0.14,0.1]), 3),
     ( 10, 0.9 , (c_double * 10)(*[0.99, 0.97, 0.95, 0.9009, 0.9, 0.61, 0.5, 0.3, 0.14,0.1]), 4),
     ( 2 , 0.9 , (c_double * 2)(*[0.61,0.5]), 0),
     ( 2 , 0.2 , (c_double * 2)(*[0.61,0.5]), 2)]
)
def test_get_num_shell_cr(p, Rs, len, expected):
    cmontecarlo_methods.get_num_shell_cr.restype = c_int64
    result = cmontecarlo_methods.get_num_shell_cr(c_double(p), Rs,len)
    assert_equal(result,expected)

@pytest.mark.parametrize(
    ['cr_idx', 'no_of_cr_shells' , 'Rs', 'expected'],
    [( 4 , 3 , (c_double * 5)(*[0.61, 0.5, 0.3, 0.14, 0.1]), 0.61),
     ( 0 , 3 , (c_double * 5)(*[0.61, 0.5, 0.3, 0.14, 0.1]), 0.61),
     ( 1 , 3 , (c_double * 5)(*[0.61, 0.5, 0.3, 0.14, 0.1]), 0.5),
     ( 3 , 3 , (c_double * 5)(*[0.61, 0.5, 0.3, 0.14, 0.1]), 0.5),
     ( 3 , 4 , (c_double * 5)(*[0.61, 0.5, 0.3, 0.14, 0.1]), 0.14),
     ( 1 , 1 , (c_double * 2)(*[0.61,0.5]), 0.61)]
)
def test_get_r(cr_idx, no_of_cr_shells, Rs, expected):
    cmontecarlo_methods.get_r.restype = c_double
    result = cmontecarlo_methods.get_r(c_int64(cr_idx), c_int64(no_of_cr_shells), Rs)
    assert_almost_equal(result, expected)
