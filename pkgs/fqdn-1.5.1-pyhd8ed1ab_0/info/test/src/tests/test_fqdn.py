# coding=utf-8
import sys

import pytest
from fqdn import FQDN


@pytest.fixture(params=(True, False))
def a_u(request):
    return request.param


class TestFQDNValidation:
    def test_constructor(self, a_u):
        with pytest.raises(ValueError):
            FQDN(None, allow_underscores=a_u)

    # Python 3-specific tests
    if sys.version_info >= (3, 0):

        def test_constructor_raises_on_bytes(self, a_u):
            with pytest.raises(ValueError):
                FQDN(b"", allow_underscores=a_u)

            with pytest.raises(ValueError):
                FQDN(b"helloworld", allow_underscores=a_u)

    def test_str(self, a_u):
        d = "greatdomain.com"
        f = FQDN(d, allow_underscores=a_u)
        assert f.absolute == str(f)

    def test_rfc_1035_s_2_3_4__label_max_length(self, a_u):
        assert FQDN(
            "www.abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijk.com",
            allow_underscores=a_u,
        ).is_valid
        assert FQDN(
            "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijk.abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijk",
            allow_underscores=a_u,
        ).is_valid

    def test_rfc_1035_s_2_3_4__label_too_long(self, a_u):
        self.__assert_invalid_from_seq("A" * 64, "com", allow_underscores=a_u)
        self.__assert_invalid_from_seq("b" * 63, "A" * 64, "com", allow_underscores=a_u)
        self.__assert_invalid_from_seq("com", "b" * 63, "A" * 64, allow_underscores=a_u)

    def test_rfc_1035_s_2_3_4__name_too_long_254_octets(self, a_u):
        parts = [(chr(ord("A") + i % 26)) for i in range(int(254 / 2) - 1)]
        parts.append("co")
        fqdn = ".".join(parts)
        assert len(fqdn) == 254
        self.__assert_invalid_from_seq(fqdn, allow_underscores=a_u)

    def test_rfc_1035_s_2_3_4__name_ok_253_octets(self, a_u):
        parts = [(chr(ord("A") + i % 26)) for i in range(int(254 / 2))]
        fqdn = ".".join(parts)
        assert len(fqdn) == 253
        self.__assert_valid_from_seq(fqdn, allow_underscores=a_u)

    def test_rfc_1035_s_3_1__trailing_byte(self, a_u):
        parts = [(chr(ord("A") + i % 26)) for i in range(int(254 / 2))]
        fqdn = ".".join(parts) + "."
        assert len(fqdn) == 254
        self.__assert_valid_from_seq(fqdn, allow_underscores=a_u)

    def test_rfc_3696_s_2__label_invalid_starts_or_ends_with_hyphen(self):
        self.__assert_invalid_fwd_and_bkwd_from_seq("-a", "com", allow_underscores=a_u)
        self.__assert_invalid_fwd_and_bkwd_from_seq("a-", "com", allow_underscores=a_u)
        self.__assert_invalid_fwd_and_bkwd_from_seq("-a-", "com", allow_underscores=a_u)

    def test_rfc_3696_s_2__preferred_form_invalid_chars(self, a_u):
        # these should use punycode instead
        self.__assert_invalid_fwd_and_bkwd_from_seq("є", "com", allow_underscores=a_u)
        self.__assert_invalid_fwd_and_bkwd_from_seq(
            "le-tour-est-joué", "com", allow_underscores=a_u
        )
        self.__assert_invalid_fwd_and_bkwd_from_seq(
            "invalid", "cóm", allow_underscores=a_u
        )
        self.__assert_invalid_fwd_and_bkwd_from_seq(
            "ich-hätte-gern-ein-Umlaut", "de", allow_underscores=a_u
        )
        self.__assert_invalid_fwd_and_bkwd_from_seq(
            "\x01", "com", allow_underscores=a_u
        )
        self.__assert_invalid_fwd_and_bkwd_from_seq(
            "x", "\x01\x02\x01", allow_underscores=a_u
        )

    def test_underscores_extra_mode(self):
        self.__assert_valid_fwd_and_bkwd_from_seq("_", "dog", allow_underscores=True)
        self.__assert_valid_fwd_and_bkwd_from_seq("i_", "dog", allow_underscores=True)
        self.__assert_valid_fwd_and_bkwd_from_seq("o_o", "dog", allow_underscores=True)

        self.__assert_invalid_fwd_and_bkwd_from_seq("_", "dog", allow_underscores=False)
        self.__assert_invalid_fwd_and_bkwd_from_seq(
            "i_", "dog", allow_underscores=False
        )
        self.__assert_invalid_fwd_and_bkwd_from_seq(
            "o_o", "dog", allow_underscores=False
        )

    def test_rfc_3696_s_2__valid(self):
        assert FQDN("net", min_labels=1, allow_underscores=a_u).is_valid
        assert FQDN("who.is", allow_underscores=a_u).is_valid
        assert FQDN("bbc.co.uk", allow_underscores=a_u).is_valid
        self.__assert_valid_fwd_and_bkwd_from_seq(
            "sh4d05-7357", "c00-mm", allow_underscores=a_u
        )

    def test_rfc_1035_s_2_3_1__label_can_have_inital_digit(self, a_u):
        self.__assert_valid_fwd_and_bkwd_from_seq("www", "1", allow_underscores=a_u)
        self.__assert_valid_fwd_and_bkwd_from_seq("1w", "1", allow_underscores=a_u)
        self.__assert_valid_fwd_and_bkwd_from_seq("1w", "a", allow_underscores=a_u)
        self.__assert_valid_fwd_and_bkwd_from_seq("1w1", "d", allow_underscores=a_u)
        self.__assert_valid_fwd_and_bkwd_from_seq("111", "a", allow_underscores=a_u)
        self.__assert_valid_fwd_and_bkwd_from_seq("www", "1a", allow_underscores=a_u)

    def test_rfc_1123__label_can_have_medial_and_terminal_digits(self, a_u):
        self.__assert_valid_fwd_and_bkwd_from_seq("www1", "a", allow_underscores=a_u)
        self.__assert_valid_fwd_and_bkwd_from_seq("ww1a", "c", allow_underscores=a_u)

        self.__assert_valid_fwd_and_bkwd_from_seq("w2w", "c", allow_underscores=a_u)
        self.__assert_valid_fwd_and_bkwd_from_seq("a111", "a", allow_underscores=a_u)
        self.__assert_valid_fwd_and_bkwd_from_seq("a1c1", "a", allow_underscores=a_u)

    def __assert_valid_fwd_and_bkwd_from_seq(self, *seq, **kwargs):
        rseq = reversed(seq)
        self.__assert_valid_from_seq(*rseq, **kwargs)

    def __assert_invalid_fwd_and_bkwd_from_seq(self, *seq, **kwargs):
        rseq = reversed(seq)
        self.__assert_invalid_from_seq(*rseq, **kwargs)

    def __assert_invalid_from_seq(self, *seq, **kwargs):
        assert not (self.__is_valid_fqdn_from_labels_seq(seq, **kwargs))

    def __assert_valid_from_seq(self, *seq, **kwargs):
        assert self.__is_valid_fqdn_from_labels_seq(seq, **kwargs)

    def __is_valid_fqdn_from_labels_seq(self, fqdn_labels_seq, **kwargs):
        fqdn = ".".join(fqdn_labels_seq)
        return FQDN(fqdn, **kwargs).is_valid


class TestMinLabels:
    def test_labels_count(self, a_u):
        assert FQDN("label").labels_count == 1
        assert FQDN("label.").labels_count == 1
        assert FQDN("label.babel").labels_count == 2
        assert FQDN("label.babel.").labels_count == 2
        assert FQDN(".label.babel.").labels_count == 3

    def test_min_labels_defaults_to_require_2(self):
        dn = FQDN("label")
        assert dn._min_labels == 2
        assert dn.labels_count == 1
        assert not dn.is_valid

    def test_min_labels_valid_set_to_1(self):
        with pytest.raises(ValueError):
            FQDN("", min_labels=1).is_valid
        assert FQDN("label", min_labels=1).is_valid
        assert not FQDN(".label", min_labels=1).is_valid
        assert FQDN("label.babel", min_labels=1).is_valid
        assert FQDN("label.babel.", min_labels=1).is_valid
        assert not FQDN(".label.babel", min_labels=1).is_valid

    def test_min_labels_valid_set_to_3(self):
        assert not FQDN("label", min_labels=3).is_valid
        assert not FQDN("label.babel", min_labels=3).is_valid
        assert not FQDN(".babel", min_labels=3).is_valid
        assert not FQDN("babel.", min_labels=3).is_valid
        assert not FQDN(".babel.", min_labels=3).is_valid
        assert not FQDN("label.babel.", min_labels=3).is_valid
        assert not FQDN(".label.babel.", min_labels=3).is_valid
        assert FQDN("fable.label.babel.", min_labels=3).is_valid
        assert FQDN("fable.label.babel", min_labels=3).is_valid


class TestAbsoluteFQDN:
    def test_absolute_fqdn(self, a_u):
        assert FQDN("trainwreck.com.", allow_underscores=a_u).is_valid_absolute is True

    def test_absolute_fqdn__fail(self, a_u):
        assert FQDN("trainwreck.com", allow_underscores=a_u).is_valid_absolute is False

    def test_to_absolute_fqdn_from_relative(self, a_u):
        assert (
            FQDN("trainwreck.com", allow_underscores=a_u).absolute == "trainwreck.com."
        )

    def test_to_absolute_fqdn_from_absolute(self, a_u):
        assert (
            FQDN("absolutetrainwreck.com.", allow_underscores=a_u).absolute
            == "absolutetrainwreck.com."
        )

    def test_to_absolute_fqdn__raises_ValueError(self, a_u):
        with pytest.raises(ValueError):
            FQDN("trainwreckcom", allow_underscores=a_u).absolute

    def test_relative_fqdn_true(self, a_u):
        assert FQDN("relative.com", allow_underscores=a_u).is_valid_relative is True

    def test_relative_fqdn_false(self, a_u):
        assert FQDN("relative.com.", allow_underscores=a_u).is_valid_relative is False


class TestRelativeFQDN:
    def test_relative_fqdn_from_relative(self, a_u):
        assert (
            FQDN("trainwreck.com", allow_underscores=a_u).relative == "trainwreck.com"
        )

    def test_relative_fqdn_from_absolute(self, a_u):
        assert (
            FQDN("trainwreck.com.", allow_underscores=a_u).relative == "trainwreck.com"
        )

    def test_relative_fqdn_from_invalid(self, a_u):
        with pytest.raises(ValueError):
            FQDN("trainwreck..", allow_underscores=a_u).relative


class TestEquality:
    def test_absolutes_are_equal(self, a_u):
        assert FQDN("trainwreck.com.", allow_underscores=a_u) == FQDN(
            "trainwreck.com.", allow_underscores=a_u
        )

    def test_relatives_are_equal(self, a_u):
        assert FQDN("trainwreck.com", allow_underscores=a_u) == FQDN(
            "trainwreck.com", allow_underscores=a_u
        )

    def test_mismatch_are_equal(self, a_u):
        assert FQDN("trainwreck.com.", allow_underscores=a_u) == FQDN(
            "trainwreck.com", allow_underscores=a_u
        )

    def test_equality_is_case_insensitive(self, a_u):
        assert FQDN(
            "all-letters-were-created-equal.com.", allow_underscores=a_u
        ) == FQDN("ALL-LETTERS-WERE-CREATED-EQUAL.COM.", allow_underscores=a_u)

    def test_strict_and_loose_can_be_equal(self):
        assert FQDN("trainwreck.com.", allow_underscores=False) == FQDN(
            "trainwreck.com", allow_underscores=True
        )


class TestHash:
    def test_is_hashable(self, a_u):
        assert hash(FQDN("trainwreck.com."))

    def test_absolutes_are_equal(self, a_u):
        assert hash(FQDN("trainwreck.com.", allow_underscores=a_u)) == hash(
            FQDN("trainwreck.com.", allow_underscores=a_u)
        )

    def test_relatives_are_equal(self, a_u):
        assert hash(FQDN("trainwreck.com", allow_underscores=a_u)) == hash(
            FQDN("trainwreck.com", allow_underscores=a_u)
        )

    def test_mismatch_are_equal(self, a_u):
        assert hash(FQDN("trainwreck.com.", allow_underscores=a_u)) == hash(
            FQDN("trainwreck.com", allow_underscores=a_u)
        )

    def test_equality_is_case_insensitive(self, a_u):
        assert hash(
            FQDN("all-letters-were-created-equal.com.", allow_underscores=a_u)
        ) == hash(FQDN("ALL-LETTERS-WERE-CREATED-EQUAL.COM.", allow_underscores=a_u))

    def test_not_equal_to_string(self, a_u):
        assert hash(FQDN("trainwreck.com.", allow_underscores=a_u)) != hash(
            "trainwreck.com."
        )

    def test_different_fqdns_are_not_equal(self, a_u):
        assert hash(FQDN("trainwreck.com.")) != hash(FQDN("test.com."))

    def test_strict_and_loose_hashs_are_equal(self):
        assert hash(FQDN("trainwreck.com.", allow_underscores=False)) == hash(
            FQDN("trainwreck.com", allow_underscores=True)
        )
