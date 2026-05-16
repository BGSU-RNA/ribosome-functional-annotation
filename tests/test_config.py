"""Unit tests for runtime configuration models (spec §7.4)."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from ribosome_state_annotator.config import CompletenessThresholds


def test_default_thresholds_match_spec_7_4() -> None:
    t = CompletenessThresholds()
    assert t.bacterial == 15
    assert t.eukaryotic == 30
    assert t.eukaryotic_organellar == 20


def test_for_classification_dispatches_correctly() -> None:
    t = CompletenessThresholds(bacterial=10, eukaryotic=20, eukaryotic_organellar=15)
    assert t.for_classification("bacterial_ribosome") == 10
    assert t.for_classification("eukaryotic_ribosome") == 20
    assert t.for_classification("eukaryotic_organellar_ribosome") == 15


def test_rejects_negative_threshold() -> None:
    with pytest.raises(ValidationError):
        CompletenessThresholds(bacterial=-1)


def test_zero_threshold_is_valid() -> None:
    """A user can disable the threshold check by setting it to 0."""
    t = CompletenessThresholds(bacterial=0)
    assert t.bacterial == 0
