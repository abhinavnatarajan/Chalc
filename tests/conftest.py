"""Shared optional dependencies for the test suite."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
	from collections.abc import Callable

	import numpy as np


@pytest.fixture(scope="session")
def bottleneck_distance() -> Callable[[np.ndarray, np.ndarray], float]:
	"""Provide GUDHI only to the tests that need it."""
	return pytest.importorskip(
		"gudhi.bottleneck",
		reason="Install the test-gudhi group to run GUDHI comparisons.",
		exc_type=ModuleNotFoundError,
	).bottleneck_distance
