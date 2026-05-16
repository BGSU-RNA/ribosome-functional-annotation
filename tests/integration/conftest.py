"""Auto-mark every test in this directory as ``@pytest.mark.network``.

The default ``pytest`` invocation excludes the ``network`` marker
(``pyproject.toml`` ``addopts = ["-m", "not network"]``). To run
the live-API tests:

::

    pytest -m network
"""

from __future__ import annotations

from collections.abc import Iterable

import pytest


def pytest_collection_modifyitems(config: pytest.Config, items: Iterable[pytest.Item]) -> None:
    """Attach the ``network`` marker to every test under ``tests/integration/``.

    Saves us from sprinkling ``pytestmark = pytest.mark.network`` at the top
    of each integration test file.
    """
    for item in items:
        if "tests/integration" in str(item.fspath):
            item.add_marker(pytest.mark.network)
