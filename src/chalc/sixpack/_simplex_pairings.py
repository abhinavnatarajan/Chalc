from __future__ import annotations

from collections.abc import Collection, Iterator
from typing import Literal, TypeVar

import numpy as np

NumRows = TypeVar("NumRows", bound=int)
NumCols = TypeVar("NumCols", bound=int)


class SimplexPairings(Collection):
	"""Persistence diagram object."""

	_paired: set[tuple[int, int]]
	_unpaired: set[int]

	@property
	def paired(self) -> set[tuple[int, int]]:
		"""Set of indices of paired simplices (read-only)."""
		return self._paired.copy()

	@property
	def unpaired(self) -> set[int]:
		"""Set of indices of unpaired simplices (read-only)."""
		return self._unpaired.copy()

	def __init__(
		self,
		paired: Collection[tuple[int, int]] = set(),
		unpaired: Collection[int] = set(),
	) -> None:
		"""Initialise the persistence diagram object."""
		self._paired = set(paired)
		self._unpaired = set(unpaired)

	def __str__(self) -> str:
		"""Return string representation of the persistence diagram."""
		return f"Paired: {self._paired}\nUnpaired: {self._unpaired}"

	def __eq__(self, other: object) -> bool:
		"""Check if two diagrams have the same paired and unpaired simplices."""
		if not isinstance(other, SimplexPairings):
			return False
		return self._paired == other._paired and self._unpaired == other._unpaired

	def __hash__(self) -> int:
		"""Return a hash of the persistence diagram."""
		return hash((self._paired, self._unpaired))

	def __bool__(self) -> bool:
		"""Return true if the diagram is non-empty."""
		# check if the diagram is non-empty
		return bool(self._paired) or bool(self._unpaired)

	# Implement the collections protocol
	def __len__(self) -> int:
		"""Return the number of features in the diagram."""
		return len(self._paired) + len(self._unpaired)

	def __iter__(self) -> Iterator[tuple[int, int] | int]:
		"""Iterate over all features in the diagram."""
		yield from self._paired
		yield from self._unpaired

	def __contains__(self, feature: object) -> bool:
		"""
		Return true if a feature is in the diagram.

		The feature to check should be either a pair of simplices (int, int)
		or a single simplex (int).
		"""
		if isinstance(feature, int) or (
			isinstance(feature, tuple) and list(map(type, feature)) == [int, int]
		):
			return feature in self._paired or feature in self._unpaired
		return False

	@classmethod
	def _from_matrices(
		cls,
		paired_matrix: np.ndarray[tuple[NumRows, Literal[2]], np.dtype[np.integer]],
		unpaired_vector: np.ndarray[tuple[int], np.dtype[np.integer]],
	) -> SimplexPairings:
		# Construct a SimplexPairings object from
		# paired and unpaired simplices represented as matrices.
		# Internal use only - used to load diagrams from file.
		paired: set[tuple[int, int]] = {
			(int(paired_matrix[i, 0]), int(paired_matrix[i, 1]))
			for i in range(paired_matrix.shape[0])
		}
		unpaired: set[int] = {int(x) for x in unpaired_vector}
		return cls(paired, unpaired)

	def paired_as_matrix(self) -> np.ndarray[tuple[int, Literal[2]], np.dtype[np.int64]]:
		"""Return a matrix representation of the finite persistence features in the diagram."""
		return np.array([[x, y] for (x, y) in self._paired], dtype=np.int64)
