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

import pickle

import pytest

from rpds import List

HASH_MSG = "Not sure List implements Hash, it has mutable methods"


def test_literalish_works():
    assert List(1, 2, 3) == List([1, 2, 3])


def test_first_and_rest():
    pl = List([1, 2])
    assert pl.first == 1
    assert pl.rest.first == 2
    assert pl.rest.rest == List()


def test_instantiate_large_list():
    assert List(range(1000)).first == 0


def test_iteration():
    assert list(List()) == []
    assert list(List([1, 2, 3])) == [1, 2, 3]


def test_push_front():
    assert List([1, 2, 3]).push_front(0) == List([0, 1, 2, 3])


def test_push_front_empty_list():
    assert List().push_front(0) == List([0])


def test_truthiness():
    assert List([1])
    assert not List()


def test_len():
    assert len(List([1, 2, 3])) == 3
    assert len(List()) == 0


def test_first_illegal_on_empty_list():
    with pytest.raises(IndexError):
        List().first


def test_rest_return_self_on_empty_list():
    assert List().rest == List()


def test_reverse():
    assert reversed(List([1, 2, 3])) == List([3, 2, 1])

    assert reversed(List()) == List()


def test_inequality():
    assert List([1, 2]) != List([1, 3])
    assert List([1, 2]) != List([1, 2, 3])
    assert List() != List([1, 2, 3])


def test_repr():
    assert str(List()) == "List([])"
    assert str(List([1, 2, 3])) in "List([1, 2, 3])"


@pytest.mark.xfail(reason=HASH_MSG)
def test_hashing():
    assert hash(List([1, 2])) == hash(List([1, 2]))
    assert hash(List([1, 2])) != hash(List([2, 1]))


def test_sequence():
    m = List("asdf")
    assert m == List(["a", "s", "d", "f"])


# Non-pyrsistent-test-suite tests


def test_drop_first():
    assert List([1, 2, 3]).drop_first() == List([2, 3])


def test_drop_first_empty():
    """
    rpds itself returns an Option<List> here but we try IndexError instead.
    """
    with pytest.raises(IndexError):
        List([]).drop_first()


def test_more_eq():
    o = object()

    assert List([o, o]) == List([o, o])
    assert List([o]) == List([o])
    assert List() == List([])
    assert not (List([1, 2]) == List([1, 3]))
    assert not (List([o]) == List([o, o]))
    assert not (List([]) == List([o]))

    assert List([1, 2]) != List([1, 3])
    assert List([o]) != List([o, o])
    assert List([]) != List([o])
    assert not (List([o, o]) != List([o, o]))
    assert not (List([o]) != List([o]))
    assert not (List() != List([]))


def test_pickle():
    assert pickle.loads(pickle.dumps(List([1, 2, 3, 4]))) == List([1, 2, 3, 4])
