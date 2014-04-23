import numpy as np
#from tardis import montecarlo
from tardis.tests import montecarlo_test_wrappers as  montecarlo
import pytest


test_line_list = np.array([10, 9, 8, 7, 6, 5, 5, 4, 3, 2, 1]).astype(np.float64)


@pytest.mark.parametrize(("insert_value", "expected_insert_position"), [
    (9.5, 0),
    (8.5, 1),
    (7.5, 2),
    (6.5, 3),
    (5.5, 4),
    (5.2, 4),
    (4.5, 6),
    (3.5, 7),
    (2.5, 8),
    (1.5, 9)])
def test_binary_search(insert_value, expected_insert_position):
    insert_position = montecarlo.binary_search_wrapper(test_line_list, insert_value, 0, len(test_line_list)-1)
    assert insert_position == expected_insert_position


@pytest.mark.parametrize(("insert_value"), [
    (10.5),
    (0.5)])
def test_binary_search_out_of_bounds(insert_value, capsys):
    with pytest.raises(ValueError):
        insert_position = montecarlo.binary_search_wrapper(test_line_list, insert_value, 0, len(test_line_list)-1)

@pytest.mark.parametrize(("insert_value", "expected_insert_position"), [
    (10.5, 0),
    (0.5, len(test_line_list))])
def test_line_search_out_of_bounds(insert_value, expected_insert_position):
    insert_position = montecarlo.line_search_wrapper(test_line_list,
                            insert_value, len(test_line_list))

    assert insert_position == expected_insert_position



    
