import pytest
from pipeline.cli.prompts import parse_stop_after, is_valid_kmer, VALID_SIZES


def test_parse_stop_after_blank_returns_none():
    assert parse_stop_after("") is None


def test_parse_stop_after_valid_int():
    assert parse_stop_after("10") == 10


def test_parse_stop_after_negative_raises():
    with pytest.raises(ValueError):
        parse_stop_after("-1")


def test_parse_stop_after_zero_raises():
    with pytest.raises(ValueError):
        parse_stop_after("0")


def test_parse_stop_after_non_int_raises():
    with pytest.raises(ValueError):
        parse_stop_after("abc")


def test_is_valid_kmer_bounds():
    assert is_valid_kmer(2) is True
    assert is_valid_kmer(6) is True
    assert is_valid_kmer(1) is False
    assert is_valid_kmer(7) is False


def test_valid_sizes_set():
    assert VALID_SIZES == {"500", "1000", "1500", "2000"}
