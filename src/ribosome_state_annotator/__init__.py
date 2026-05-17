"""Annotate functional states of complete bacterial and eukaryotic ribosome assemblies.

Public API surface — re-exported from :mod:`.api`, :mod:`.models`, and
:mod:`.exceptions` so callers can ``from ribosome_state_annotator import …``
without reaching into submodules.
"""

from __future__ import annotations

from ribosome_state_annotator.api import (
    annotate_assembly,
    annotate_many,
    annotate_pdb,
)
from ribosome_state_annotator.config import CompletenessThresholds
from ribosome_state_annotator.exceptions import (
    ApiRequestError,
    CoordinateDownloadError,
    CoordinateParsingError,
    CorrespondenceMappingError,
    RibosomeAnnotatorError,
    UnsupportedRibosomeError,
)
from ribosome_state_annotator.models import (
    Anticodon,
    AnticodonResidue,
    AssemblyContext,
    BasePair,
    ChainRef,
    Codon,
    CodonResidue,
    CorrespondenceResult,
    LargeScaleMovements,
    LigandRef,
    RibosomeAnnotation,
    RibosomeClassification,
    RibosomeStatus,
    TRNAmRNAInteraction,
)

__version__ = "0.1.0.dev0"

__all__ = [
    "Anticodon",
    "AnticodonResidue",
    "ApiRequestError",
    "AssemblyContext",
    "BasePair",
    "ChainRef",
    "Codon",
    "CodonResidue",
    "CompletenessThresholds",
    "CoordinateDownloadError",
    "CoordinateParsingError",
    "CorrespondenceMappingError",
    "CorrespondenceResult",
    "LargeScaleMovements",
    "LigandRef",
    "RibosomeAnnotation",
    "RibosomeAnnotatorError",
    "RibosomeClassification",
    "RibosomeStatus",
    "TRNAmRNAInteraction",
    "UnsupportedRibosomeError",
    "__version__",
    "annotate_assembly",
    "annotate_many",
    "annotate_pdb",
]
