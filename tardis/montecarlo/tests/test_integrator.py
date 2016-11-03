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
import numpy as np

from ctypes import CDLL, byref, c_uint, c_int, c_int64, c_double, c_ulong, POINTER
from numpy.testing import assert_equal, assert_almost_equal, assert_approx_equal

from tardis import __path__ as path
from tardis.montecarlo.struct import IndexPair

cmontecarlo_filepath = os.path.join(path[0], 'montecarlo', 'montecarlo.so')
cmontecarlo_methods = CDLL(cmontecarlo_filepath)

RS_SIX_SHELL = np.array([2.24640000e+15, 2.07792000e+15, 1.90944000e+15, 1.74096000e+15, 1.57248000e+15, 1.40400000e+15, 1.23552000e+15])


class DataCollection:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


@pytest.fixture
def six_shell_three_lines():
    rs = np.array([2.24640000e+15, 2.07792000e+15, 1.90944000e+15, 1.74096000e+15, 1.57248000e+15, 1.40400000e+15, 1.23552000e+15])
    line_nu = np.array([ 6.16510725e+14, 4.56675108e+14, 1.59835930e+14])
    s_ul = np.array([[  1.03081070e-05,   6.88799038e-06,   5.34986499e-06, 4.08517213e-06,   3.23265549e-06,   2.69619598e-06],
                     [  1.61772204e-05,   1.06697121e-05,   8.03779370e-06, 6.42442807e-06,   5.19508084e-06,   4.39577502e-06],
                     [  8.32617749e-06,   5.70817366e-06,   4.52953255e-06, 4.26894544e-06,   3.28933528e-06,   2.75121093e-06]])
    taus = np.array([[  2.62870433,   2.62870433,   2.62870433,   2.62870433, 2.62870433,   2.62870433],
                     [ 19.81254555,  19.81254555,  19.81254555,  19.81254555, 19.81254555,  19.81254555],
                     [  3.11753516,   3.11753516,   3.11753516,   3.11753516, 3.11753516,   3.11753516]])
    return DataCollection(Rs    = (c_double * rs.size)(*rs.flatten(order='C')),
                          lines = (c_double * line_nu.size)(*line_nu.flatten(order='C')),
                          s_ul  = (c_double * s_ul.size)(*s_ul.flatten(order='C')),
                          taus  = (c_double * taus.size)(*taus.flatten(order='C')),
                          Rmax  = rs.max(),
                          llen  = line_nu.size)
@pytest.fixture
def six_shell_one_lines():
    rs = np.array([2.24640000e+15, 2.07792000e+15, 1.90944000e+15, 1.74096000e+15, 1.57248000e+15, 1.40400000e+15, 1.23552000e+15])
    line_nu = np.array([  4.566751e+14 ])
    s_ul = np.array([1.61772204e-05,   1.06697121e-05,   8.03779370e-06, 6.42442807e-06,   5.19508084e-06,   4.39577502e-06])
    taus = np.array([19.812546, 19.812546, 19.812546, 19.812546, 19.812546, 19.812546])
    return DataCollection(Rs    = (c_double * rs.size)(*rs.flatten(order='C')),
                          lines = (c_double * line_nu.size)(*line_nu.flatten(order='C')),
                          s_ul  = (c_double * s_ul.size)(*s_ul.flatten(order='C')),
                          taus  = (c_double * taus.size)(*taus.flatten(order='C')),
                          Rmax  = rs.max(),
                          llen  = line_nu.size)

@pytest.fixture
def six_shell_four_lines():
    rs = np.array([2.24640000e+15, 2.07792000e+15, 1.90944000e+15, 1.74096000e+15, 1.57248000e+15, 1.40400000e+15, 1.23552000e+15])
    line_nu = np.array([4.56684917e+14, 4.56684013e+14, 4.56680743e+14, 4.56675108e+14])
    s_ul = np.array([ [  1.61906752e-05,   1.06801494e-05,   8.02856211e-06, 6.42414065e-06,   5.21304216e-06,   4.36988692e-06],
                    [  1.60534857e-05,   1.05708064e-05,   7.95511056e-06, 6.35116471e-06,   5.15967121e-06,   4.32514816e-06],
                    [  1.45385557e-05,   9.59885565e-06,   7.20776792e-06, 5.79331823e-06,   4.66891071e-06,   3.92746368e-06],
                    [  1.61905785e-05,   1.06693121e-05,   8.03749634e-06, 6.41695541e-06,   5.20709427e-06,   4.36995683e-06]])
    taus = np.array([[11.020823, 11.020823, 11.020823, 11.020823, 11.020823, 11.020823],
                     [ 4.580072,  4.580072,  4.580072,  4.580072,  4.580072,  4.580072],
                     [ 2.290083,  2.290083,  2.290083,  2.290083,  2.290083,  2.290083],
                     [19.812545, 19.812545, 19.812545, 19.812545, 19.812545, 19.812545]])

    return DataCollection(Rs    = (c_double * rs.size)(*rs.flatten(order='C')),
                          lines = (c_double * line_nu.size)(*line_nu.flatten(order='C')),
                          s_ul  = (c_double * s_ul.size)(*s_ul.flatten(order='C')),
                          taus  = (c_double * taus.size)(*taus.flatten(order='C')),
                          Rmax  = rs.max(),                        
                          llen  = line_nu.size)
                          

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
    ['no_of_shells', 'p' , 'R_ph', 'expected'],
    [( 3 , 1.0 , 2.0, 4),
     ( 1 , 1.0 , 2.0, 2),
     ( 1 , 2.0 , 1.0, 0),
     ( 3 , 1797120000000000.0 , 1235520000000000.0, 0),
     ( 3 ,  898560000000000.0 , 1235520000000000.0, 4),
    ])
def test_get_cr_start(no_of_shells, p, R_ph, expected):
    cmontecarlo_methods.get_cr_start.restype = c_int64
    result = cmontecarlo_methods.get_cr_start(c_int64(no_of_shells), c_double(p), c_double(R_ph))
    assert_equal(result,expected)

@pytest.mark.parametrize(
    ['cr_idx', 'no_of_cr_shells', 'expected'],
    [( 0  ,  1  ,  0),
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
     ( 5 , 0.09, (c_double * 5)(*[0.61, 0.5, 0.3, 0.14,0.1]), 5),
     ( 10, 0.9 , (c_double * 10)(*[0.99, 0.97, 0.95, 0.9009, 0.9, 0.61, 0.5, 0.3, 0.14,0.1]), 4),
     ( 2 , 0.9 , (c_double * 2)(*[0.61,0.5]), 0),
     ( 2 , 0.2 , (c_double * 2)(*[0.61,0.5]), 2)
     ])
def test_get_num_shell_cr(p, Rs, len, expected):
    cmontecarlo_methods.get_num_shell_cr.restype = c_int64
    result = cmontecarlo_methods.get_num_shell_cr(c_double(p), Rs,len)
    assert_equal(result,expected)

@pytest.mark.parametrize(
    ['cr_idx', 'no_of_cr_shells' , 'Rs', 'expected'],
    [( 4 , 3 , (c_double * 5)(*[0.61, 0.5, 0.3, 0.14, 0.1]), 0.5),
     ( 5 , 3 , (c_double * 5)(*[0.61, 0.5, 0.3, 0.14, 0.1]), 0.61),
     ( 0 , 3 , (c_double * 5)(*[0.61, 0.5, 0.3, 0.14, 0.1]), 0.61),
     ( 1 , 3 , (c_double * 5)(*[0.61, 0.5, 0.3, 0.14, 0.1]), 0.5),
     ( 2 , 3 , (c_double * 5)(*[0.61, 0.5, 0.3, 0.14, 0.1]), 0.3),
     ( 3 , 3 , (c_double * 5)(*[0.61, 0.5, 0.3, 0.14, 0.1]), 0.3),
     ( 3 , 4 , (c_double * 5)(*[0.61, 0.5, 0.3, 0.14, 0.1]), 0.14),
     ( 1 , 1 , (c_double * 2)(*[0.61,0.5]), 0.61)]
)
def test_get_r(cr_idx, no_of_cr_shells, Rs, expected):
    cmontecarlo_methods.get_r.restype = c_double
    result = cmontecarlo_methods.get_r(c_int64(cr_idx), c_int64(no_of_cr_shells), Rs)
    assert_almost_equal(result, expected)

@pytest.mark.parametrize(
    ['nu','rmax_frac','cr_idx', 'no_of_cr_shells', 'inv_ct', 'input_data', 'expected'],
    [(c_double( 6.16510725e+14),0.71, 0, 4, c_double(2.969765804826852e-17),
        six_shell_three_lines, (645473954328893.38,640896021318693.62) ),
    (c_double(6.16510725e+14), 0.71, 3, 4, c_double(2.969765804826852e-17),
        six_shell_three_lines, (629288929493385.5, 603732520506614.38) ),
#    (c_double( 4.56675108e+14), c_double(RS_SIX_SHELL.max()*0.88), 3, 4, c_double(2.969765804826852e-17),
#        (c_double * 7)(*RS_SIX_SHELL), #Rs 
#        (c_double * 4)(*[4.56684917e+14, 4.56684013e+14, 4.56680743e+14, 4.56675108e+14]), 4, # Line_nu, linelen
#        (0,3)),    
    ]
)
def test_test_nu_limits_for_crossing_and_p(nu, rmax_frac, cr_idx, no_of_cr_shells, inv_ct, input_data, expected):
    indata = input_data()
    p = c_double(indata.Rmax*rmax_frac)
    cmontecarlo_methods.test_nu_limits_for_crossing_and_p.restype = c_double
    blu = cmontecarlo_methods.test_nu_limits_for_crossing_and_p(nu, p, cr_idx, no_of_cr_shells, inv_ct,
            indata.Rs, indata.lines, indata.llen, 0)
    red = cmontecarlo_methods.test_nu_limits_for_crossing_and_p(nu, p, cr_idx, no_of_cr_shells, inv_ct, 
            indata.Rs, indata.lines, indata.llen, 1)
    assert_approx_equal(blu, expected[0],significant=9)
    assert_approx_equal(red, expected[1],significant=9)

@pytest.mark.parametrize(
    ['nu_blu', 'nu_red' , 'input_data', 'expected'],
    [( c_double(5.0e+14), c_double(4.0e+14), six_shell_four_lines, (0,3) ),
     ( c_double(5.0e+14), c_double(4.56684013e+14), six_shell_four_lines, (0,1) ),
     ( c_double(7.0e+14), c_double(4.566751e+14), six_shell_one_lines, (0,0) ),
     ( c_double(7.0e+14), c_double(4.566751e+14), six_shell_one_lines, (0,0) ),
     ( c_double(4.59596043e+14), c_double(4.44756413e+14), six_shell_three_lines, (1,1) ),
])
def test_nu_idx_from_nu_pair(capfd,nu_blu, nu_red, input_data, expected):
    indata = input_data()
    cmontecarlo_methods.nu_idx_from_nu_pair.restype = IndexPair
    result = cmontecarlo_methods.nu_idx_from_nu_pair(nu_blu, nu_red, indata.lines, indata.llen)

    out = capfd.readouterr()
    print out[0]

    assert_equal(result.start, expected[0])
    assert_equal(result.end,   expected[1])

@pytest.mark.parametrize(
    ['nu','p','cr_idx', 'no_of_cr_shells', 'inv_ct',  'input_data', 'expected'],
    [(c_double(6.16510725e+14), 0.71, 3, 4, c_double(2.969765804826852e-17),
        six_shell_three_lines, (0,0) ),
     (c_double(4.56675108e+14), 0.71, 3, 4, c_double(2.969765804826852e-17),
        six_shell_three_lines, (1,1) ),
     (c_double(4.56680743e+14), 0.71, 3, 4, c_double(2.969765804826852e-17),
        six_shell_four_lines, (0,3) ),
    ])
def test_find_nu_limits_for_crossing_and_p(capfd, nu, p, cr_idx, no_of_cr_shells, inv_ct, input_data, expected):
    indata = input_data()
    p = c_double(indata.Rmax*p)
    cmontecarlo_methods.find_nu_limits_for_crossing_and_p.restype = IndexPair
    result = cmontecarlo_methods.find_nu_limits_for_crossing_and_p(nu, p, cr_idx, no_of_cr_shells, inv_ct, 
        indata.Rs, indata.lines, indata.llen)

    out = capfd.readouterr()
    print out[0]
    
#    assert_equal(result.start, expected[0])
#    assert_equal(result.end,   expected[1])

@pytest.mark.parametrize(
    ['I_nu','nu_lims','sh_idx','len','expected'],
    [(0.0, IndexPair(0,0), 0, 3,  1.03081070e-05),
     (0.0, IndexPair(0,0), 1, 3,  4.08517213e-06),
#     (0.0, IndexPair(0,1), 1, 3,  4.08517213e-06)
])
def test_sum_lines(I_nu,nu_lims, sh_idx, len, expected, atom_data=six_shell_three_lines()):
    cmontecarlo_methods.sum_lines.restype = c_double
    result = cmontecarlo_methods.sum_lines(nu_lims, c_double(I_nu), atom_data.taus, atom_data.s_ul, sh_idx, len)
#    assert_almost_equal(result,expected)

@pytest.mark.parametrize(
    ['len','Rmax'],
    [(20, 2.24640000e+15),]
   )

def test_integrate_intensity(len,Rmax):
    ps = np.linspace(0.999, 0, num = len)*Rmax
    I_nu = np.array([2.40739908129e+27, 2.40039294821e+27, 2.78750654685e+27, 3.03077239309e+27, 3.02260473791e+27, 3.30099043765e+27, 3.29209678908e+27, 3.28285455902e+27, 3.434666461e+27, 3.2633625314e+27, 3.25313169884e+27, 2.46201222899e+27, 1.87332590148e+27, 1.00422409649e+27, 7.53168073762e+26, 5.37977197133e+26, 3.58651466602e+26, 2.05262250637e+27, 2.22838279668e+27, 2.21568602442e+27])
    I_p = I_nu*ps
    expected = 8 * np.pi**2 *  np.trapz(y = I_p[::-1],x = ps[::-1])
    ps   = (c_double * len)(*ps)
    I_nu = (c_double * len)(*I_nu)
    cmontecarlo_methods.integrate_intensity.restype = c_double
    result = cmontecarlo_methods.integrate_intensity(I_nu,ps,len)
    assert_approx_equal(result,expected,significant=9)



@pytest.mark.parametrize(
    ['i','j','len','expected'],
    [(0, 0, 6, 1.03081070e-05),
     (1, 1, 6, 1.06697121e-05),
     (5, 0, 6, 2.69619598e-06),
     (3, 2, 6, 4.26894544e-06),
    ]
)
def test_intex_macro(i,j,len,expected,atom_data=six_shell_three_lines()):
    cmontecarlo_methods.test_index_macro.restype = c_double
    result = cmontecarlo_methods.test_index_macro(atom_data.s_ul,i,j,len)
    assert_almost_equal(result,expected)

