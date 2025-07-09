"""Implementation of a 6-pack of persistence diagrams."""

from __future__ import annotations

from collections.abc import (
	Callable,
	Collection,
	ItemsView,
	Iterator,
	KeysView,
	Mapping,
	Sequence,
	ValuesView,
)
from typing import Literal, get_args, overload

import numpy as np
from h5py import Dataset, Group

type DiagramName = Literal["ker", "cok", "dom", "cod", "im", "rel"]
"""Names of diagrams in a 6-pack of persistence diagrams."""


class SixPack(Mapping):
	"""6-pack of persistence diagrams."""

	__slots__ = ("_dimensions", "_entrance_times", "_simplex_pairings")
	_dimensions: np.ndarray[tuple[int], np.dtype[np.int64]]
	_entrance_times: np.ndarray[tuple[int], np.dtype[np.float64]]
	_simplex_pairings: dict[DiagramName, SimplexPairings]

	@classmethod
	def names(cls) -> tuple[DiagramName, ...]:
		"""Return the names of the diagrams in the 6-pack."""
		return get_args(DiagramName.__value__)

	@property
	def entrance_times(self) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
		"""Entrance times of the simplices."""
		temp = self._entrance_times[:]
		temp.setflags(write=False)
		return temp

	@property
	def dimensions(self) -> np.ndarray[tuple[int], np.dtype[np.int64]]:
		"""Dimensions of the simplices."""
		temp = self._dimensions[:]
		temp.setflags(write=False)
		return temp

	def max_nonempty_dimension(self) -> int:
		"""Get the maximum dimension of features across all diagrams.

		Returns -1 if there are no features in the 6-pack.
		"""
		max_dim = -1
		for name, dgm in self.items():
			if not dgm:
				continue
			if name != "ker":
				max_dim = max(
					max_dim,
					*(self._dimensions[simplex] for simplex in dgm.unpaired),
					*(self._dimensions[s1] for s1, _ in dgm.paired),
				)
			else:
				max_dim = max(
					max_dim,
					*(self._dimensions[simplex] - 1 for simplex in dgm.unpaired),
					*(self._dimensions[s1] - 1 for s1, _ in dgm.paired),
				)
		return max_dim

	def __init__(
		self,
		kernel: SimplexPairings | None = None,
		cokernel: SimplexPairings | None = None,
		domain: SimplexPairings | None = None,
		codomain: SimplexPairings | None = None,
		image: SimplexPairings | None = None,
		relative: SimplexPairings | None = None,
		entrance_times: Sequence[float] = [],
		dimensions: Sequence[int] = [],
	) -> None:
		"""Initialise the 6-pack of persistence diagrams."""
		self._simplex_pairings = {}

		self._simplex_pairings["ker"] = kernel or SimplexPairings()
		self._simplex_pairings["cok"] = cokernel or SimplexPairings()
		self._simplex_pairings["dom"] = domain or SimplexPairings()
		self._simplex_pairings["cod"] = codomain or SimplexPairings()
		self._simplex_pairings["im"] = image or SimplexPairings()
		self._simplex_pairings["rel"] = relative or SimplexPairings()
		self._entrance_times = np.array(
			entrance_times,
		)
		self._dimensions = np.array(dimensions)

	def __getitem__(self, key: DiagramName) -> SimplexPairings:
		"""Access a specific diagram in the 6-pack."""
		return self._simplex_pairings[key]

	def items(self) -> ItemsView[DiagramName, SimplexPairings]:
		"""View of the diagrams in the 6-pack."""
		return self._simplex_pairings.items()

	def __contains__(self, key: object) -> bool:
		"""Return true if a diagram is in the 6-pack."""
		return key in self._simplex_pairings

	def __iter__(
		self,
	) -> Iterator[DiagramName]:
		"""Iterate over all diagrams in the 6-pack."""
		yield from self._simplex_pairings

	def keys(self) -> KeysView[DiagramName]:
		"""View of the names of the diagrams in the 6-pack."""
		return self._simplex_pairings.keys()

	def values(self) -> ValuesView[SimplexPairings]:
		"""View of the diagrams in the 6-pack."""
		return self._simplex_pairings.values()

	def __len__(self) -> int:
		"""Return the number of diagrams in the 6-pack."""
		return 6

	def __bool__(self) -> bool:
		"""Return true if any diagram in the 6-pack is non-empty."""
		return any(self._simplex_pairings.values())

	def __eq__(self, other: object) -> bool:
		"""Check if two 6-packs of persistence diagrams are identical."""
		if not isinstance(other, SixPack):
			return False
		return (
			all(self[name] == other[name] for name in self)
			and all(self._entrance_times == other._entrance_times)
			and all(self._dimensions == other._dimensions)
		)

	def __hash__(self) -> int:
		"""Return a hash of the 6-pack of persistence diagrams."""
		return hash(
			(
				self._simplex_pairings.items(),
				self._entrance_times,
				self._dimensions,
			),
		)

	def num_features(self) -> int:
		"""Count the total number of features across all diagrams in the 6-pack."""
		return sum(len(dgm) for dgm in self.values())

	def filter(self, func: Callable[[str, int, float, float], bool]) -> SixPack:
		"""Filter out features in the diagram.

		``func`` should take four arguments: the name of the diagram to which a feature
		belongs, the dimension of the feature, and its birth and death times, and should return
		a boolean indicating whether the feature should be kept.
		"""
		for name, diagram in self.items():
			paired = frozenset(
				(birth_simplex, death_simplex)
				for birth_simplex, death_simplex in diagram.paired
				if func(
					name,
					self._dimensions[birth_simplex] - (1 if name == "ker" else 0),
					self._entrance_times[birth_simplex],
					self._entrance_times[death_simplex],
				)
			)
			unpaired = frozenset(
				simplex
				for simplex in diagram.unpaired
				if func(
					name,
					self._dimensions[simplex] - (1 if name == "ker" else 0),
					self._entrance_times[simplex],
					np.inf,
				)
			)
			diagram._paired = paired  # noqa: SLF001
			diagram._unpaired = unpaired  # noqa: SLF001
		self._cleanup()
		return self

	def _cleanup(self) -> SixPack:
		"""Remove all simplices not associated with a feature."""
		# Find all simplices that are part of some feature
		active_simplices: set[int] = set()
		for pairings in self.values():
			for simplex_1, simplex_2 in pairings.paired:
				active_simplices.add(simplex_1)
				active_simplices.add(simplex_2)
			for simplex in pairings.unpaired:
				active_simplices.add(simplex)

		# Create a map from old indices to new indices
		sorted_active_simplices = sorted(active_simplices)
		index_map = {old: new for new, old in enumerate(sorted_active_simplices)}

		# Filter dimensions and entrance times
		if active_simplices:
			self._dimensions = self._dimensions[sorted_active_simplices]
			self._entrance_times = self._entrance_times[sorted_active_simplices]
		else:
			self._dimensions = np.array([], dtype=np.int64)
			self._entrance_times = np.array([], dtype=np.float64)

		# Re-index the simplex pairings
		for name in self._simplex_pairings:
			pairings = self._simplex_pairings[name]
			new_paired = frozenset((index_map[s1], index_map[s2]) for s1, s2 in pairings.paired)
			new_unpaired = frozenset(index_map[s] for s in pairings.unpaired)
			self._simplex_pairings[name] = SimplexPairings(new_paired, new_unpaired)

		return self

	@overload
	def get_matrix(
		self,
		diagram_name: DiagramName,
		dimension: int,
	) -> np.ndarray[tuple[int, Literal[2]], np.dtype[np.float64]]: ...

	@overload
	def get_matrix(
		self,
		diagram_name: DiagramName,
		dimension: Sequence[int] | None = None,
	) -> list[np.ndarray[tuple[int, Literal[2]], np.dtype[np.float64]]]: ...

	def get_matrix(self, diagram_name, dimension=None):
		r"""
		Get a specific diagram as a matrix of birth and death times.

		Args:
			diagram_name :
				One of ``'ker'``, ``'cok'``, ``'dom'``,
				``'cod'``, ``'im'``, or ``'rel'``.
			dimension    :
				Dimension(s) of the diagram desired.
				If a list is provided then a list of matrices is returned,
				with the order of matrices respecting the order of entries of `dim`.
				If `dimension` is not provided then the returned matrix will contain
				persistent features from all homological dimensions
				from zero to ``max(self.dimensions)``.

		Returns:
			An :math:`m \times 2` matrix whose rows are
			a pair of birth and death times, or a list of such matrices.

		"""
		if dimension is None:
			simplices = set()
			dgm = self[diagram_name]
			for b, d in dgm.paired:
				simplices.add(b)
				simplices.add(d)
			for s in dgm.unpaired:
				simplices.add(s)
			relevant_dims = np.array(self._dimensions)[list(simplices)]
			if relevant_dims.size == 0 or relevant_dims.ndim == 0:
				# there are no features in the diagram
				return []
			dimension = list(range(max([x.item() for x in relevant_dims])))

		if isinstance(dimension, int):
			# a p-dimensional homology class is captured by
			# a pairing of (p+1) simplices for kernels
			if diagram_name == "ker":
				dimension += 1
			diagram = self[diagram_name]
			entrance_times = np.array(self._entrance_times)
			paired_list = np.array(
				[pair for pair in diagram.paired if self._dimensions[pair[0]] == dimension],
				dtype=np.int64,
			)
			unpaired_list = np.array(
				[sigma for sigma in diagram.unpaired if self._dimensions[sigma] == dimension],
				dtype=np.int64,
			)
			paired_matrix = entrance_times[paired_list.flatten()].reshape((-1, 2))
			unpaired_matrix = entrance_times[unpaired_list][:, np.newaxis]
			inf_array = np.array([np.inf for _ in range(len(unpaired_list))])[:, np.newaxis]
			unpaired_matrix = np.concatenate((unpaired_matrix, inf_array), axis=1)
			return np.concatenate((paired_matrix, unpaired_matrix), axis=0)

		if isinstance(dimension, Sequence):
			return [self.get_matrix(diagram_name, d) for d in dimension]

		errmsg = "dim must be an integer or a list of integers"
		raise TypeError(errmsg)

	def save(self, file: Group) -> None:
		"""
		Save a 6-pack of persistence diagrams to a HDF5 file or group.

		Args:
			file : A h5py file or group.

		"""
		file["entrance_times"] = np.array(self._entrance_times, dtype=np.float64)
		file["dimensions"] = np.array(self._dimensions, dtype=np.int64)
		for dgm_name in self:
			grp = file.require_group(dgm_name)
			grp["unpaired"] = np.array(list(self[dgm_name].unpaired), dtype=np.int64)
			grp["paired"] = self[dgm_name]._paired_as_matrix()

	@classmethod
	def from_file(cls, file: Group) -> SixPack:
		"""
		Load a 6-pack of persistence diagrams from a HDF5 file or group.

		Args:
			file : A h5py file or group.

		"""
		dgms = SixPack()
		if any(
			attr not in file
			for attr in ["entrance_times", "dimensions", "ker", "cok", "dom", "cod", "im", "rel"]
		):
			errmsg = "Invalid file: missing fields"
			raise RuntimeError(errmsg)

		dgms._entrance_times = (
			entrance_times[:]
			if isinstance((entrance_times := file["entrance_times"]), Dataset)
			and entrance_times[:].dtype == np.float64
			else np.array([], dtype=np.float64)
		)
		dgms._dimensions = (
			dimensions[:]
			if isinstance((dimensions := file["dimensions"]), Dataset)
			and dimensions[:].dtype == np.int64
			else np.array([], dtype=np.int64)
		)
		names: tuple[DiagramName, ...] = ("ker", "cok", "dom", "cod", "im", "rel")
		for diagram_name in names:
			grp = file[diagram_name]
			if isinstance(grp, Group):
				temp = {}
				for name in ("paired", "unpaired"):
					if name not in grp:
						errmsg = f"Invalid file: missing field {diagram_name}.{name}."
						raise RuntimeError(errmsg)
					if not isinstance(grp[name], Dataset):
						errmsg = f"Invalid file: expected {diagram_name}.{name} "
						f"to be a Dataset(dtype=np.int64), but found {type(grp[name])}."
						raise TypeError(errmsg)
					temp[name] = grp[name]
					if temp[name].dtype != np.int64:
						errmsg = f"Invalid file: expected {diagram_name}.{name} "
						"to be a Dataset(dtype=np.int64), "
						f"but found Dataset(dtype={temp[name].dtype})."
						raise RuntimeError(errmsg)
				paired = temp["paired"].astype(int)
				unpaired = temp["unpaired"].astype(int)
				dgms._simplex_pairings[diagram_name] = SimplexPairings._from_matrices(  # noqa: SLF001
					paired[...],
					unpaired[...],
				)
			else:
				errmsg = f"Invalid file: expected {diagram_name} to be a group "
				f"but found a {type(grp)}."
				raise TypeError(errmsg)

		return dgms


class SimplexPairings(Collection):
	"""Persistence diagram represented by paired and unpaired simplices."""

	__slots__ = ("_paired", "_unpaired")
	_paired: frozenset[tuple[int, int]]
	_unpaired: frozenset[int]

	@property
	def paired(self) -> frozenset[tuple[int, int]]:
		"""Set of indices of paired simplices."""
		return self._paired

	@property
	def unpaired(self) -> frozenset[int]:
		"""Set of indices of unpaired simplices."""
		return self._unpaired

	def __init__(
		self,
		paired: Collection[tuple[int, int]] = frozenset(),
		unpaired: Collection[int] = frozenset(),
	) -> None:
		"""Initialise a diagram with its paired and unpaired simplices."""
		self._paired = frozenset(paired)
		self._unpaired = frozenset(unpaired)

	def __str__(self) -> str:
		"""Represent the persistence diagram as a string."""
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
		if isinstance(feature, int):
			return feature in self._unpaired
		if isinstance(feature, tuple) and tuple(map(type, feature)) == (int, int):
			return feature in self._paired
		return False

	@classmethod
	def _from_matrices[NRows: int](
		cls,
		paired_matrix: np.ndarray[tuple[NRows, Literal[2]], np.dtype[np.integer]],
		unpaired_vector: np.ndarray[tuple[int], np.dtype[np.integer]],
	) -> SimplexPairings:
		# Construct a SimplexPairings object from
		# paired and unpaired simplices represented as matrices.
		# Internal use only - used to load diagrams from file.
		paired = frozenset(
			tuple(x.item() for x in paired_matrix[i, :]) for i in range(paired_matrix.shape[0])
		)
		unpaired = frozenset(x.item() for x in unpaired_vector)
		return cls(paired, unpaired)

	def _paired_as_matrix(self) -> np.ndarray[tuple[int, Literal[2]], np.dtype[np.int64]]:
		r"""Return the pairings of simplex indices in a :math:`N\times 2` matrix."""
		return np.array([[x, y] for (x, y) in self._paired], dtype=np.int64)
