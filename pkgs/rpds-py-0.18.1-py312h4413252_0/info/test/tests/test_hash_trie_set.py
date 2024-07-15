"""
Modified from the pyrsistent test suite.

Pre-modification, these were MIT licensed, and are copyright:

    Copyright (c) 2022 Tobias Gustafsson

    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use,
    copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following
    conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
    OTHER DEALINGS IN THE SOFTWARE.
"""

from collections import abc
import pickle

import pytest

from rpds import HashTrieSet

HASH_MSG = "Not sure HashTrieSet implements Hash, it has mutable methods"


def test_key_is_tuple():
    with pytest.raises(KeyError):
        HashTrieSet().remove((1, 1))


def test_key_is_not_tuple():
    with pytest.raises(KeyError):
        HashTrieSet().remove("asdf")


@pytest.mark.xfail(reason=HASH_MSG)
def test_supports_hash():
    assert hash(HashTrieSet((1, 2))) == hash(HashTrieSet(1, 2))


def test_empty_truthiness():
    assert HashTrieSet([1])
    assert not HashTrieSet()


def test_contains_elements_that_it_was_initialized_with():
    initial = [1, 2, 3]
    s = HashTrieSet(initial)

    assert set(s) == set(initial)
    assert len(s) == len(set(initial))


def test_is_immutable():
    s1 = HashTrieSet([1])
    s2 = s1.insert(2)

    assert s1 == HashTrieSet([1])
    assert s2 == HashTrieSet([1, 2])

    s3 = s2.remove(1)
    assert s2 == HashTrieSet([1, 2])
    assert s3 == HashTrieSet([2])


def test_remove_when_not_present():
    s1 = HashTrieSet([1, 2, 3])
    with pytest.raises(KeyError):
        s1.remove(4)


def test_discard():
    s1 = HashTrieSet((1, 2, 3))
    assert s1.discard(3) == HashTrieSet((1, 2))
    assert s1.discard(4) == s1


def test_is_iterable():
    assert sum(HashTrieSet([1, 2, 3])) == 6


def test_contains():
    s = HashTrieSet([1, 2, 3])

    assert 2 in s
    assert 4 not in s


def test_supports_set_operations():
    s1 = HashTrieSet([1, 2, 3])
    s2 = HashTrieSet([3, 4, 5])

    assert s1 | s2 == HashTrieSet([1, 2, 3, 4, 5])
    assert s1.union(s2) == s1 | s2

    assert s1 & s2 == HashTrieSet([3])
    assert s1.intersection(s2) == s1 & s2

    assert s1 - s2 == HashTrieSet([1, 2])
    assert s1.difference(s2) == s1 - s2

    assert s1 ^ s2 == HashTrieSet([1, 2, 4, 5])
    assert s1.symmetric_difference(s2) == s1 ^ s2


def test_supports_set_comparisons():
    s1 = HashTrieSet([1, 2, 3])
    s3 = HashTrieSet([1, 2])
    s4 = HashTrieSet([1, 2, 3])

    assert HashTrieSet([1, 2, 3, 3, 5]) == HashTrieSet([1, 2, 3, 5])
    assert s1 != s3

    assert s3 < s1
    assert s3 <= s1
    assert s3 <= s4

    assert s1 > s3
    assert s1 >= s3
    assert s4 >= s3


def test_repr():
    rep = repr(HashTrieSet([1, 2]))
    assert rep == "HashTrieSet({1, 2})" or rep == "HashTrieSet({2, 1})"

    rep = repr(HashTrieSet(["1", "2"]))
    assert rep == "HashTrieSet({'1', '2'})" or rep == "HashTrieSet({'2', '1'})"


def test_update():
    assert HashTrieSet([1, 2, 3]).update([3, 4, 4, 5]) == HashTrieSet(
        [1, 2, 3, 4, 5],
    )


def test_update_no_elements():
    s1 = HashTrieSet([1, 2])
    assert s1.update([]) == s1


def test_iterable():
    assert HashTrieSet(iter("a")) == HashTrieSet(iter("a"))


def test_more_eq():
    # Non-pyrsistent-test-suite test
    o = object()

    assert HashTrieSet([o]) == HashTrieSet([o])
    assert HashTrieSet([o, o]) == HashTrieSet([o, o])
    assert HashTrieSet([o]) == HashTrieSet([o, o])
    assert HashTrieSet() == HashTrieSet([])
    assert not (HashTrieSet([1, 2]) == HashTrieSet([1, 3]))
    assert not (HashTrieSet([o, 1]) == HashTrieSet([o, o]))
    assert not (HashTrieSet([]) == HashTrieSet([o]))

    assert HashTrieSet([1, 2]) != HashTrieSet([1, 3])
    assert HashTrieSet([]) != HashTrieSet([o])
    assert not (HashTrieSet([o]) != HashTrieSet([o]))
    assert not (HashTrieSet([o, o]) != HashTrieSet([o, o]))
    assert not (HashTrieSet([o]) != HashTrieSet([o, o]))
    assert not (HashTrieSet() != HashTrieSet([]))

    assert HashTrieSet([1, 2]) == {1, 2}
    assert HashTrieSet([1, 2]) != {1, 2, 3}
    assert HashTrieSet([1, 2]) != [1, 2]


def test_more_set_comparisons():
    s = HashTrieSet([1, 2, 3])

    assert s == s
    assert not (s < s)
    assert s <= s
    assert not (s > s)
    assert s >= s


def test_pickle():
    assert pickle.loads(
        pickle.dumps(HashTrieSet([1, 2, 3, 4])),
    ) == HashTrieSet([1, 2, 3, 4])


def test_instance_of_set():
    assert isinstance(HashTrieSet(), abc.Set)


def test_lt_le_gt_ge():
    assert HashTrieSet({}) < {1}
    assert HashTrieSet({}) <= {1}
    assert HashTrieSet({1}) > set()
    assert HashTrieSet({1}) >= set()
