from tardis import util


def test_atomic_number2element_symbol():
    assert util.atomic_number2element_symbol(14) == 'Si'

def test_element_symbol2atomic_number():
    assert util.element_symbol2atomic_number('Si') == 14