"""Routines for computing 6-packs of persistence diagrams."""

from __future__ import annotations

from collections.abc import Collection, ItemsView, KeysView, Mapping, Sequence, ValuesView
from typing import Callable, Iterator, Literal, TypeVar, overload

import numpy as np
from h5py import Dataset, Group
from phimaker import compute_ensemble

from chalc.chromatic import alpha, delcech, delrips
from chalc.filtration import FilteredComplex

__all__ = [
	"from_filtration",
	"compute",
	"SimplexPairings",
	"DiagramEnsemble",
]

ChromaticMethod = {
	"alpha": alpha,
	"delcech": delcech,
	"delrips": delrips,
}

BdMatType = list[tuple[bool, int, list[int]]]

UnknownType = TypeVar("UnknownType")


# bitmask that represents a list of colours
def _colours_to_bitmask(colours: Collection[int]) -> int:
	return sum(2**i for i in colours)


# list of colours represented by a bitmask
def _bitmask_to_colours(b: int) -> list[int]:
	i = 0
	res = []
	while b:
		if b & 1:
			res.append(i)
		i += 1
		b >>= 1
	return res


# number of colours specified by a bitmask
def _num_colours_in_bitmask(b: int) -> int:
	return len(_bitmask_to_colours(b))


# check if the colours specified by bitmask b1 are a subset of those specified by b2
def _colours_are_subset(b1: int, b2: int) -> bool:
	return not (~b2 & b1)


def _get_diagrams(
	K: FilteredComplex,
	check_in_domain: Callable[[int], bool],
	max_dgm_dim: int,
	tolerance: float = 0,
) -> DiagramEnsemble:
	# Build up the matrix
	matrix: BdMatType = []
	entrance_times: list[float] = []
	dimensions: list[int] = []
	for column in K.serialised():
		facet_idxs = column[0]
		dimension = max(0, len(facet_idxs) - 1)
		# Only need k-skeleton for k <= max_dgm_dim + 1
		if dimension > max_dgm_dim + 1:
			break
		dimensions.append(dimension)
		entrance_times.append(column[2])
		matrix.append((check_in_domain(column[3]), dimension, facet_idxs))
	d = compute_ensemble(matrix)
	dgms = DiagramEnsemble(
		SimplexPairings(d.ker.paired, d.ker.unpaired),
		SimplexPairings(d.cok.paired, d.cok.unpaired),
		SimplexPairings(d.g.paired, d.g.unpaired),
		SimplexPairings(d.f.paired, d.f.unpaired),
		SimplexPairings(d.im.paired, d.im.unpaired),
		SimplexPairings(d.rel.paired, d.rel.unpaired),
		entrance_times,
		dimensions,
	)
	return dgms.threshold(tolerance)


def from_filtration(
	K: FilteredComplex,
	dom: Collection[int] | int | None = None,
	k: int | None = None,
	max_diagram_dimension: int | None = None,
	tolerance: float = 0,
) -> DiagramEnsemble:
	"""
	Compute the 6-pack of persistence diagrams associated to a filtered simplicial complex with coloured vertices.

	Given a filtered chromatic simplicial complex :math:`K` and a subcomplex :math:`L` of :math:`K`, this function computes the 6-pack of persistence diagram associated with the inclusion map :math:`f : L \\hookrightarrow K`. The subcomplex is specified by
	the colours of its vertices, or by an integer :math:`k` wherein all simplices with :math:`k` or fewer colours are considered part of the subcomplex.

	Args:
		K                     : A filtered chromatic simplicial complex.
		dom                   : Integer or collection of integers describing the colours of the points in the domain (the subcomplex :math:`L`).
		k                     : If not ``None``, then the domain is taken to be the :math:`k`-chromatic subcomplex of :math:`K`, i.e., the subcomplex of simplices having at most :math:`k` colours.
		max_diagram_dimension : Maximum homological dimension for which the persistence diagrams are computed. By default diagrams of all dimensions are computed.
		tolerance             : Retain only points with persistence strictly greater than this value.

	Returns:
		Diagrams corresponding to the following persistence modules (where :math:`H_*` is the persistent homology functor and :math:`f_*` is the induced map on persistent homology):

		#. :math:`H_*(L)` (domain)
		#. :math:`H_*(K)` (codomain)
		#. :math:`\\ker(f_*)` (kernel)
		#. :math:`\\mathrm{coker}(f_*)` (cokernel)
		#. :math:`\\mathrm{im}(f_*)` (image)
		#. :math:`H_*(K, L)` (relative homology)

		Each diagram is represented by sets of paired and unpaired simplices,
		and contain simplices of all dimensions. ``dgms`` also contains the
		entrance times of the simplices and their dimensions.
	"""
	if k is None and dom is None:
		raise RuntimeError("At least one of k or dom must be provided")
	if k is not None and dom is not None:
		raise RuntimeError("Only one of k or dom is allowed")
	if dom is not None:
		if isinstance(dom, int):
			dom = [
				dom,
			]
		colours_bitmask = _colours_to_bitmask(dom)

		def check_in_domain(b):
			return _colours_are_subset(b, colours_bitmask)
	elif k is not None:

		def check_in_domain(b):
			return _num_colours_in_bitmask(b) <= k  # k-chromatic simplex

	# if(K, L) is a pair of simplicial complexes where dim(L) <= dim(K) = d
	# then all of the sixpack diagrams except for the relative homology
	# can be non-zero upto dimension d, and relative can be nontrivial upto dim = d+1
	if max_diagram_dimension is None:
		max_diagram_dimension = K.dimension + 1
	else:
		max_diagram_dimension = min(max_diagram_dimension, K.dimension + 1)
	return _get_diagrams(K, check_in_domain, max_diagram_dimension, tolerance)


def compute(
	x: np.ndarray[tuple[int, int], np.dtype[np.float64]],
	colours: Sequence[int],
	dom: Collection[int] | int | None = None,
	k: int | None = None,
	method: Literal["alpha", "delcech", "delrips"] = "alpha",
	max_diagram_dimension: int | None = None,
	tolerance: float = 0,
) -> DiagramEnsemble:
	"""
	Compute the 6-pack of persistence diagrams of a coloured point-cloud.

	This function constructs a filtered simplicial complex :math:`K` from the point cloud, and computes the 6-pack of persistence diagrams associated with the inclusion :math:`f : L \\hookrightarrow K` where :math:`L` is some filtered subcomplex of :math:`K`.

	Args:
		x                     : Numpy matrix whose columns are points.
		colours               : Sequence of integers describing the colours of the points.
		dom                   : Integer or collection of integers describing the colours of the points in the domain (the subcomplex :math:`L`).
		k                     : If not ``None``, then the domain is taken to be the :math:`k`-chromatic subcomplex of :math:`K`, i.e., the subcomplex of simplices having at most :math:`k` colours.
		method                : Filtration used to construct the chromatic complex. Must be one of ``'alpha'``, ``'delcech'``, or ``'delrips'``.
		max_diagram_dimension : Maximum homological dimension for which the persistence diagrams are computed. By default diagrams of all dimensions are computed.
		tolerance             : Retain only points with persistence strictly greater than this value.

	Returns :
		Diagrams corresponding to the following persistence modules (where :math:`H_*` is the persistent homology functor and :math:`f_*` is the induced map on persistent homology):

		#. :math:`H_*(L)` (domain)
		#. :math:`H_*(K)` (codomain)
		#. :math:`\\ker(f_*)` (kernel)
		#. :math:`\\mathrm{coker}(f_*)` (cokernel)
		#. :math:`\\mathrm{im}(f_*)` (image)
		#. :math:`H_*(K, L)` (relative homology)

		Each diagram is represented by sets of paired and unpaired simplices,
		and contains simplices of all dimensions. ``dgms`` also contains the
		entrance times of the simplices and their dimensions.
	"""
	if k is None and dom is None:
		raise RuntimeError("At least one of k or dom must be provided")
	if k is not None and dom is not None:
		raise RuntimeError("Only one of k or dom is allowed")
	if dom is not None:
		# new colours: 1 -> domain, 0 -> codomain
		if isinstance(dom, int):
			dom = [
				dom,
			]
		new_colours = list(np.isin(colours, list(dom)).astype(np.int64))

		def check_in_domain(b: int):
			return b == 2
	elif k is not None:

		def check_in_domain(b: int):
			return _num_colours_in_bitmask(b) <= k  # k-chromatic simplex

		new_colours = list(colours)
	# can have non-zero homology upto dimension d for everything except rel
	# rel can have non-zero homology in dimension d+1
	if max_diagram_dimension is None:
		max_diagram_dimension = x.shape[0]
	else:
		max_diagram_dimension = min(max_diagram_dimension, x.shape[0])
	# Compute chromatic complex
	if method in ChromaticMethod.keys():
		K, _ = ChromaticMethod[method](x, new_colours)
	else:
		raise RuntimeError(f"method must be one of {ChromaticMethod.keys()}")
	return _get_diagrams(K, check_in_domain, max_diagram_dimension, tolerance)


class SimplexPairings(Collection):
	"""
	Persistence diagram object, represented by a
	set of simplex pairings and a set of unpaired simplices.
	"""

	_paired: set[tuple[int, int]]
	_unpaired: set[int]

	@property
	def paired(self) -> set[tuple[int, int]]:
		"""
		Set of indices of paired simplices (read-only).
		"""
		return self._paired.copy()

	@property
	def unpaired(self) -> set[int]:
		"""
		Set of indices of unpaired simplices (read-only).
		"""
		return self._unpaired.copy()

	def __init__(
		self, paired: Collection[tuple[int, int]] = set(), unpaired: Collection[int] = set()
	) -> None:
		self._paired = set(paired)
		self._unpaired = set(unpaired)

	def __str__(self):
		return f"Paired: {self._paired}\nUnpaired: {self._unpaired}"

	def __eq__(self, other: SimplexPairings) -> bool:
		return self._paired == other._paired and self._unpaired == other._unpaired

	def __bool__(self) -> bool:
		"""
		Returns true if the diagram is non-empty.
		"""
		# check if the diagram is non-empty
		return bool(self._paired) or bool(self._unpaired)

	# Implement the collections protocol
	def __len__(self) -> int:
		"""
		Returns the number of features in the diagram.
		"""
		# return the number of features in the diagram
		return len(self._paired) + len(self._unpaired)

	def __iter__(self) -> Iterator[tuple[int, int] | int]:
		"""
		Iterate over all features in the diagram.
		"""
		yield from self._paired
		yield from self._unpaired

	def __contains__(self, feature: tuple[int, int] | int) -> bool:
		"""
		Returns true if a feature, represented by either a pair of simplices or a single unpaired simplex, is in the diagram.
		"""
		return feature in self._paired or feature in self._unpaired

	@classmethod
	def _from_matrices(
		cls,
		paired_matrix: np.ndarray[tuple[int, Literal[2]], np.dtype[np.int64]],
		unpaired_vector: np.ndarray[int, np.dtype[np.int64]],
	) -> SimplexPairings:
		"""
		Construct a SimplexPairings object from paired and unpaired simplices represented as matrices.
		Internal use only - used to load diagrams from file.
		"""
		paired: set[tuple[int, int]] = {
			(int(paired_matrix[i, 0]), int(paired_matrix[i, 1]))
			for i in range(paired_matrix.shape[0])
		}
		unpaired: set[int] = set(int(x) for x in unpaired_vector)
		return cls(paired, unpaired)

	def paired_as_matrix(self) -> np.ndarray[tuple[int, Literal[2]], np.dtype[np.int64]]:
		"""
		Returns a matrix representation of the finite persistence features in the diagram.
		"""
		return np.concatenate(list(self._paired)).reshape((-1, 2))


class DiagramEnsemble(Mapping):
	"""
	6-pack of persistence diagrams.
	"""

	DiagramName = Literal["ker", "cok", "dom", "cod", "im", "rel"]
	"""Type hint for the names of the diagrams in the 6-pack."""
	diagram_names: list[DiagramName] = ["ker", "cok", "dom", "cod", "im", "rel"]
	"""Names of the diagrams in the 6-pack."""

	@property
	def entrance_times(self) -> np.ndarray[int, np.dtype[np.float64]]:
		"""
		Entrance times of the simplices.
		"""
		temp = self._entrance_times[:]
		temp.setflags(write=False)
		return temp

	@property
	def dimensions(self) -> np.ndarray[int, np.dtype[np.int64]]:
		"""
		Dimensions of the simplices.
		"""
		temp = self._dimensions[:]
		temp.setflags(write=False)
		return temp

	def __init__(
		self,
		ker: SimplexPairings = SimplexPairings(),
		cok: SimplexPairings = SimplexPairings(),
		dom: SimplexPairings = SimplexPairings(),
		cod: SimplexPairings = SimplexPairings(),
		im: SimplexPairings = SimplexPairings(),
		rel: SimplexPairings = SimplexPairings(),
		entrance_times: Sequence[float] = [],
		dimensions: Sequence[int] = [],
	):
		self._simplex_pairings: dict[DiagramEnsemble.DiagramName, SimplexPairings] = dict()

		self._simplex_pairings["ker"] = ker
		self._simplex_pairings["cok"] = cok
		self._simplex_pairings["dom"] = dom
		self._simplex_pairings["cod"] = cod
		self._simplex_pairings["im"] = im
		self._simplex_pairings["rel"] = rel
		self._entrance_times: np.ndarray[int, np.dtype[np.float64]] = np.array(entrance_times)
		self._dimensions: np.ndarray[int, np.dtype[np.int64]] = np.array(dimensions)

	def __getitem__(self, key: DiagramName) -> SimplexPairings:
		"""
		Access a specific diagram in the 6-pack.
		"""
		return self._simplex_pairings[key]

	def get(self, key: DiagramName, default: UnknownType = None) -> SimplexPairings | UnknownType:
		"""
		Access a specific diagram in the 6-pack, with a default value if the diagram does not exist.
		"""
		return self._simplex_pairings.get(key, default)

	def items(self) -> ItemsView[DiagramName, SimplexPairings]:
		"""
		View of the diagrams in the 6-pack.
		"""
		return self._simplex_pairings.items()

	def __contains__(self, key: DiagramName) -> bool:
		return key in self._simplex_pairings.keys()

	def __iter__(
		self,
	) -> Iterator[tuple[DiagramName, SimplexPairings]]:
		"""
		Iterate over all diagrams in the 6-pack.
		"""
		yield from zip(self._simplex_pairings.keys(), self._simplex_pairings.values())

	def keys(self) -> KeysView[DiagramName]:
		"""
		View of the names of the diagrams in the 6-pack.
		"""
		return self._simplex_pairings.keys()

	def values(self) -> ValuesView[SimplexPairings]:
		"""
		View of the diagrams in the 6-pack.
		"""
		return self._simplex_pairings.values()

	def __len__(self) -> int:
		return 6

	def __bool__(self) -> bool:
		"""
		Returns true if any diagram in the 6-pack is non-empty.
		"""
		return any(self._simplex_pairings.values())

	def __eq__(self, other: DiagramEnsemble) -> bool:
		return (
			all(self[name] == other[name] for name in DiagramEnsemble.diagram_names)
			and all(self._entrance_times == other._entrance_times)
			and all(self._dimensions == other._dimensions)
		)

	def num_features(self) -> int:
		"""
		Count the total number of features across all diagrams in the 6-pack.
		"""
		return sum(len(dgm) for dgm in self.values())

	def threshold(self, tolerance: float) -> DiagramEnsemble:
		"""
		Discard all features with persistence ``<=tolerance``.
		"""
		for diagram_name in self.keys():
			pairings = {
				(b, d)
				for b, d in self[diagram_name]._paired
				if self._entrance_times[d] - self._entrance_times[b] > tolerance
			}
			self._simplex_pairings[diagram_name]._paired = pairings
		return self

	@overload
	def get_matrix(
		self, diagram_name: DiagramName, dim: int
	) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]: ...

	@overload
	def get_matrix(
		self,
		diagram_name: DiagramName,
		dim: list[int] | None = None,
	) -> list[np.ndarray[tuple[int, int], np.dtype[np.float64]]]: ...

	def get_matrix(self, diagram_name, dim=None):
		"""
		Get a specific diagram as a matrix of birth and death times.

		Args:
			diagram_name : One of ``'ker'``, ``'cok'``, ``'dom'``, ``'cod'``, ``'im'``, or ``'rel'``.
			dim          : Dimension(s) of the diagram desired. If a list is provided then a list of matrices is returned, with the order of matrices respecting the order of entries of `dim`. If `dim` is not provided then the returned matrix will contain persistent features from all homological dimensions from zero to ``max(self.dimensions)``.

		Returns:
			An :math:`m \\times 2` matrix whose rows are a pair of birth and death times, or a list of such matrices.
		"""
		if dim is None:
			simplices = set()
			dgm = self[diagram_name]
			for b, d in dgm._paired:
				simplices.add(b)
				simplices.add(d)
			for s in dgm._unpaired:
				simplices.add(s)
			dim = list(range(max(np.array(self._dimensions)[list(simplices)])))
		if isinstance(dim, int):
			# a p-dimensional homology class is captured by a pairing of (p+1) simplices for kernels
			if diagram_name == "ker":
				dim += 1
			diagram = self[diagram_name]
			entrance_times = np.array(self._entrance_times)
			paired_list = np.array(
				[pair for pair in diagram._paired if self._dimensions[pair[0]] == dim],
				dtype=np.int64,
			)
			unpaired_list = np.array(
				[sigma for sigma in diagram._unpaired if self._dimensions[sigma] == dim],
				dtype=np.int64,
			)
			paired_matrix = entrance_times[paired_list.flatten()].reshape((-1, 2))
			unpaired_matrix = entrance_times[unpaired_list][:, np.newaxis]
			inf_array = np.array([np.inf] * len(unpaired_list))[:, np.newaxis]
			unpaired_matrix = np.concatenate((unpaired_matrix, inf_array), axis=1)
			return np.concatenate((paired_matrix, unpaired_matrix), axis=0)
		elif isinstance(dim, list):
			return [self.get_matrix(diagram_name, d) for d in dim]

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
			grp["unpaired"] = np.array(list(diagram._unpaired), dtype=np.int64)
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
			raise RuntimeError("Invalid file: missing fields")

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
				temp = dict()
				for name in ["paired", "unpaired"]:
					if name not in grp:
						raise RuntimeError(f"Invalid file: missing field {diagram_name}.{name}.")
					if not isinstance(grp[name], Dataset):
						raise RuntimeError(
							f"Invalid file: expected {diagram_name}.{name} to be a Dataset(dtype=np.int64), but found {type(grp[name])}."
						)
					temp[name] = grp[name]
					if temp[name].dtype != np.int64:
						raise RuntimeError(
							f"Invalid file: expected {diagram_name}.{name} to be a Dataset(dtype=np.int64), but found Dataset(dtype={temp[name].dtype})."
						)
				paired = temp["paired"].astype(int)
				unpaired = temp["unpaired"].astype(int)
				dgms._simplex_pairings[diagram_name] = SimplexPairings._from_matrices(
					paired[...], unpaired[...]
				)
			else:
				raise RuntimeError(
					f"Invalid file: expected {diagram_name} to be a group but found a {type(grp)}."
				)

		return dgms
