"""Typed exception hierarchy (spec §16).

Most "this assembly isn't supported" cases are represented as structured
``RibosomeAnnotation(status="skipped", skip_reason=...)`` results rather
than exceptions. The exceptions here are reserved for **unexpected**
failures: a remote API returns an unparseable payload, a coordinate file
can't be parsed, the package itself is invoked with invalid arguments.
"""

from __future__ import annotations


class RibosomeAnnotatorError(Exception):
    """Base class for all errors raised by this package."""


class ApiRequestError(RibosomeAnnotatorError):
    """Raised when an external HTTP call to RCSB or BGSU fails or returns garbage."""


class CoordinateDownloadError(RibosomeAnnotatorError):
    """Raised when a biological-assembly mmCIF file cannot be retrieved."""


class CoordinateParsingError(RibosomeAnnotatorError):
    """Raised when Gemmi cannot parse a coordinate file."""


class CorrespondenceMappingError(RibosomeAnnotatorError):
    """Raised when the BGSU correspondence response is malformed or empty."""


class UnsupportedRibosomeError(RibosomeAnnotatorError):
    """Raised for hard-failure unsupported cases.

    Most unsupported assemblies (NMR, archaeal, partial, etc.) are returned
    as ``status="skipped"`` annotation results. This exception is reserved
    for the rare cases where a structured skip is not possible, e.g. an
    invalid PDB ID format passed directly to the library API.
    """
