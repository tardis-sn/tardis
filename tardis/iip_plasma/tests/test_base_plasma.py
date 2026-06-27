from __future__ import annotations

from tardis.iip_plasma.base import BasePlasma
from tardis.iip_plasma.properties.base import Input, ProcessingPlasmaProperty


class InputA(Input):
    outputs = ("a",)


class InputB(Input):
    outputs = ("b",)


class SumAB(ProcessingPlasmaProperty):
    outputs = ("sum_ab",)

    def calculate(self, a: float, b: float) -> float:
        return a + b


class TwiceSum(ProcessingPlasmaProperty):
    outputs = ("twice_sum",)

    def calculate(self, sum_ab: float) -> float:
        return 2.0 * sum_ab


def test_resolve_update_list_uses_cached_independent_copy() -> None:
    plasma = BasePlasma([InputA, InputB, SumAB, TwiceSum], a=1.0, b=2.0)

    update_list = plasma._resolve_update_list(["a", "b"])

    assert update_list == ["SumAB", "TwiceSum"]
    assert plasma._update_list_cache[frozenset(("a", "b"))] == update_list

    update_list.append("mutated")

    assert plasma._resolve_update_list(["b", "a"]) == ["SumAB", "TwiceSum"]
