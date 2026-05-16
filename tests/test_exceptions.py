"""Unit tests for the exception hierarchy (spec §16)."""

from __future__ import annotations

import pytest

from ribosome_state_annotator.exceptions import (
    ApiRequestError,
    CoordinateDownloadError,
    CoordinateParsingError,
    CorrespondenceMappingError,
    RibosomeAnnotatorError,
    UnsupportedRibosomeError,
)

ALL_SUBCLASSES = [
    ApiRequestError,
    CoordinateDownloadError,
    CoordinateParsingError,
    CorrespondenceMappingError,
    UnsupportedRibosomeError,
]


def test_base_inherits_from_exception() -> None:
    assert issubclass(RibosomeAnnotatorError, Exception)


@pytest.mark.parametrize("cls", ALL_SUBCLASSES)
def test_subclasses_inherit_from_base(cls: type[Exception]) -> None:
    assert issubclass(cls, RibosomeAnnotatorError)


@pytest.mark.parametrize("cls", ALL_SUBCLASSES)
def test_subclasses_can_be_raised_and_caught_as_base(cls: type[Exception]) -> None:
    with pytest.raises(RibosomeAnnotatorError):
        raise cls("synthetic failure")


def test_distinct_subclasses_are_not_each_other() -> None:
    """A specific subclass should not be caught by an unrelated sibling subclass."""
    with pytest.raises(ApiRequestError):
        try:
            raise ApiRequestError("network down")
        except CoordinateParsingError:
            pytest.fail("ApiRequestError should not be caught as CoordinateParsingError")
