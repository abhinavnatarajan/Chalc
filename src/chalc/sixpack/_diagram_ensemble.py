from __future__ import annotations

from collections.abc import (
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

from ._simplex_pairings import SimplexPairings


class DiagramEnsemble(Mapping):
	"""6-pack of persistence diagrams."""

	type DiagramName = Literal["ker", "cok", "dom", "cod", "im", "rel"]
	"""Names of the diagrams in the 6-pack."""

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

	def __init__(
		self,
		ker: SimplexPairings | None = None,
		cok: SimplexPairings | None = None,
		dom: SimplexPairings | None = None,
		cod: SimplexPairings | None = None,
		im: SimplexPairings | None = None,
		rel: SimplexPairings | None = None,
		entrance_times: Sequence[float] = [],
		dimensions: Sequence[int] = [],
	) -> None:
		"""Initialise the 6-pack of persistence diagrams."""
		self._simplex_pairings: dict[DiagramEnsemble.DiagramName, SimplexPairings] = {}

		self._simplex_pairings["ker"] = ker or SimplexPairings()
		self._simplex_pairings["cok"] = cok or SimplexPairings()
		self._simplex_pairings["dom"] = dom or SimplexPairings()
		self._simplex_pairings["cod"] = cod or SimplexPairings()
		self._simplex_pairings["im"] = im or SimplexPairings()
		self._simplex_pairings["rel"] = rel or SimplexPairings()
		self._entrance_times: np.ndarray[tuple[int], np.dtype[np.float64]] = np.array(
			entrance_times,
		)
		self._dimensions: np.ndarray[tuple[int], np.dtype[np.int64]] = np.array(dimensions)

	def __getitem__(self, key: DiagramName) -> SimplexPairings:
		"""Access a specific diagram in the 6-pack."""
		return self._simplex_pairings[key]

	def get[U](self, key: DiagramName, default: U = None) -> SimplexPairings | U:
		"""
		Access a specific diagram in the 6-pack.

		Returns a default value if the diagram does not exist.
		"""
		return self._simplex_pairings.get(key, default)

	def items(self) -> ItemsView[DiagramName, SimplexPairings]:
		"""View of the diagrams in the 6-pack."""
		return self._simplex_pairings.items()

	def __contains__(self, key: object) -> bool:
		"""Return true if a diagram is in the 6-pack."""
		if not isinstance(key, str):
			return False
		return key in self._simplex_pairings

	def __iter__(
		self,
	) -> Iterator[tuple[DiagramName, SimplexPairings]]:
		"""Iterate over all diagrams in the 6-pack."""
		yield from zip(self._simplex_pairings.keys(), self._simplex_pairings.values(), strict=True)

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
		"""Check if two 6-packs of persistence diagrams are equal."""
		if not isinstance(other, DiagramEnsemble):
			return False
		return (
			all(
				self[name] == other[name]
				for name in get_args(DiagramEnsemble.DiagramName.__value__)
			)
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

	def threshold(self, tolerance: float) -> DiagramEnsemble:
		"""Discard all features with persistence ``<=tolerance``."""
		for diagram_name in self.keys():
			pairings = {
				(b, d)
				for b, d in self[diagram_name].paired
				if self._entrance_times[d] - self._entrance_times[b] > tolerance
			}
			self._simplex_pairings[diagram_name]._paired = pairings  # noqa: SLF001
		return self

	@overload
	def get_matrix(
		self,
		diagram_name: DiagramName,
		dim: int,
	) -> np.ndarray[tuple[int, Literal[2]], np.dtype[np.float64]]: ...

	@overload
	def get_matrix(
		self,
		diagram_name: DiagramName,
		dim: list[int] | None = None,
	) -> list[np.ndarray[tuple[int, Literal[2]], np.dtype[np.float64]]]: ...

	def get_matrix(self, diagram_name, dim=None):
		r"""
		Get a specific diagram as a matrix of birth and death times.

		Args:
			diagram_name :
				One of ``'ker'``, ``'cok'``, ``'dom'``,
				``'cod'``, ``'im'``, or ``'rel'``.
			dim          :
				Dimension(s) of the diagram desired.
				If a list is provided then a list of matrices is returned,
				with the order of matrices respecting the order of entries of `dim`.
				If `dim` is not provided then the returned matrix will contain
				persistent features from all homological dimensions
				from zero to ``max(self.dimensions)``.

		Returns:
			An :math:`m \times 2` matrix whose rows are
			a pair of birth and death times, or a list of such matrices.

		"""
		if dim is None:
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
			dim = list(range(max([x.item() for x in relevant_dims])))

		if isinstance(dim, int):
			# a p-dimensional homology class is captured by
			# a pairing of (p+1) simplices for kernels
			if diagram_name == "ker":
				dim += 1
			diagram = self[diagram_name]
			entrance_times = np.array(self._entrance_times)
			paired_list = np.array(
				[pair for pair in diagram.paired if self._dimensions[pair[0]] == dim],
				dtype=np.int64,
			)
			unpaired_list = np.array(
				[sigma for sigma in diagram.unpaired if self._dimensions[sigma] == dim],
				dtype=np.int64,
			)
			paired_matrix = entrance_times[paired_list.flatten()].reshape((-1, 2))
			unpaired_matrix = entrance_times[unpaired_list][:, np.newaxis]
			inf_array = np.array([np.inf for _ in range(len(unpaired_list))])[:, np.newaxis]
			unpaired_matrix = np.concatenate((unpaired_matrix, inf_array), axis=1)
			return np.concatenate((paired_matrix, unpaired_matrix), axis=0)

		if isinstance(dim, list):
			return [self.get_matrix(diagram_name, d) for d in dim]

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
		for name, diagram in self:
			grp = file.require_group(name)
			grp["unpaired"] = np.array(list(diagram.unpaired), dtype=np.int64)
			grp["paired"] = diagram.paired_as_matrix()

	@classmethod
	def from_file(cls, file: Group) -> DiagramEnsemble:
		"""
		Load a 6-pack of persistence diagrams from a HDF5 file or group.

		Args:
			file : A h5py file or group.

		"""
		dgms = DiagramEnsemble()
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
		names: list[DiagramEnsemble.DiagramName] = ["ker", "cok", "dom", "cod", "im", "rel"]
		for diagram_name in names:
			grp = file[diagram_name]
			if isinstance(grp, Group):
				temp = {}
				for name in ["paired", "unpaired"]:
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
