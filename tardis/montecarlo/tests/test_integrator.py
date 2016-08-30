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
from ctypes import CDLL, byref, c_uint, c_int64, c_double, c_ulong, POINTER
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
    result = cmontecarlo_methods.get_cr_sign(cr_idx, no_of_cr_shells)
    assert_equal(result,expected)

@pytest.mark.parametrize(
    ['no_of_cr_shells', 'p' , 'R_ph', 'expected'],
    [( 3 , 1797120000000000.0 , 1235520000000000.0, 0),
     ( 3 ,  898560000000000.0 , 1235520000000000.0, 3),
     ( 3 , 1235520000000000.0 , 1235520000000000.0, 0),
)
def test_get_cr_start(no_of_cr_shells, p, R_ph, expected):
    result = cmontecarlo_methods.get_cr_start(no_of_cr_shells,p,R_ph)
    assert_equal(result,expected)

