"""Runtime configuration models (spec §7.4).

The package keeps biological reference data in :mod:`constants` (constants
that should never change at runtime) and tunable knobs here (values the
user may override per-invocation through the CLI or library API).
"""

from __future__ import annotations

import logging

from pydantic import BaseModel, Field

from ribosome_state_annotator.models import RibosomeClassification

logger = logging.getLogger(__name__)


class CompletenessThresholds(BaseModel):
    """Minimum ribosomal-protein chain counts per classification (spec §7.4).

    A real complete ribosome carries 50+ protein chains; these thresholds
    are deliberately conservative — they exist to exclude obviously
    incomplete depositions, not to police biological completeness. Users
    can raise (or lower) them per-invocation if they have stricter or more
    permissive needs.
    """

    bacterial: int = Field(default=15, ge=0)
    eukaryotic: int = Field(default=30, ge=0)
    eukaryotic_organellar: int = Field(default=20, ge=0)

    def for_classification(self, classification: RibosomeClassification) -> int:
        """Return the threshold that applies to ``classification``."""
        if classification == "bacterial_ribosome":
            return self.bacterial
        if classification == "eukaryotic_ribosome":
            return self.eukaryotic
        # eukaryotic_organellar_ribosome — exhaustive over the Literal.
        return self.eukaryotic_organellar
