"""Routines for computing 6-packs of persistence diagrams."""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from collections.abc import (
	Collection,
	ItemsView,
	Iterator,
	KeysView,
	Mapping,
	Sequence,
	ValuesView,
)
from itertools import combinations
from typing import TYPE_CHECKING, Literal, NewType, Self, get_args, overload

import numpy as np
from h5py import Dataset, Group
from phimaker import compute_ensemble, compute_ensemble_cylinder

from chalc.chromatic import alpha, delcech, delrips

if TYPE_CHECKING:
	from collections.abc import (
		Sequence,
	)
	from typing import TypeIs

	from chalc.filtration import Filtration


def is_homogenous_type[T](x: Collection, typ: type[T]) -> TypeIs[Collection[T]]:
	return all(isinstance(i, typ) for i in x)


__all__ = [
	"FiltrationInclusion",
	"FiltrationMorphism",
	"KChromaticInclusion",
	"KChromaticQuotient",
	"SimplexPairings",
	"SixPack",
	"SubChromaticInclusion",
	"compute",
	"from_filtration",
]


class FiltrationMorphism(ABC):
	"""Represents a map of filtrations.

	This is an abstract class. To specify a morphism, instantiate
	one of the subclasses.
	"""

	def compute_diagrams(
		self,
		filtration: Filtration,
		max_dgm_dim: int,
	) -> SixPack:
		"""Compute the 6-pack of persistence diagrams induced by this morphism."""
		return self._compute_diagrams(filtration, max_dgm_dim)

	@abstractmethod
	def _compute_diagrams(
		self,
		filtration: Filtration,
		max_dgm_dim: int,
	) -> SixPack:
		pass


class FiltrationInclusion(FiltrationMorphism, ABC):
	"""Represents the inclusion of a filtered subcomplex."""

	def simplex_in_domain(
		self,
		column: tuple[list[int], int, float, list[int]],
		filtration: Filtration,
	) -> bool:
		"""Check if a simplex is in the domain of the inclusion map.

		A simplex is identified by its column in the boundary matrix of the filtration.
		"""
		return self._simplex_in_domain(column, filtration)

	@abstractmethod
	def _simplex_in_domain(
		self,
		column: tuple[list[int], int, float, list[int]],
		filtration: Filtration,
	) -> bool:
		pass

	def _compute_diagrams(
		self,
		filtration: Filtration,
		max_dgm_dim: int,
	) -> SixPack:
		# Build up the matrix
		matrix: list[tuple[bool, int, list[int]]] = []  # [(in_domain, dimension, facet_idxs)]
		entrance_times: list[float] = []
		dimensions: list[int] = []
		for column in filtration.serialised():
			facet_idxs = column[0]
			dimension = max(0, len(facet_idxs) - 1)
			# Only need k-skeleton for k <= max_dgm_dim + 1
			if dimension > max_dgm_dim + 1:
				break
			entrance_time = column[2]
			matrix.append(
				(
					self._simplex_in_domain(column, filtration),
					dimension,
					facet_idxs,
				),
			)
			dimensions.append(dimension)
			entrance_times.append(entrance_time)

		logger = logging.getLogger(__name__)
		logger.debug("Boundary matrix: \n%s", matrix)

		d = compute_ensemble(matrix)
		return SixPack(
			SimplexPairings(d.ker.paired, d.ker.unpaired),
			SimplexPairings(d.cok.paired, d.cok.unpaired),
			SimplexPairings(d.g.paired, d.g.unpaired),
			SimplexPairings(d.f.paired, d.f.unpaired),
			SimplexPairings(d.im.paired, d.im.unpaired),
			SimplexPairings(d.rel.paired, d.rel.unpaired),
			entrance_times,
			dimensions,
		)


class SubChromaticInclusion(tuple, FiltrationInclusion):
	"""Corresponds to the inclusion of a particular subset or subsets of colours.

	If the subset is a collection of integers, then all simplices whose colours
	lie in that subset are included. Otherwise, if each element of the subset
	is itself a collection of integers, then all simplices whose colours lie in
	any of those collections are included.

	For example, ``SubChromaticInclusion(combinations(all_colours, k))`` represents the same
	morphism as ``KChromaticInclusion(k)``.
	"""

	__slots__ = ()

	def __new__(cls, x: Collection[Collection[int]] | Collection[int]) -> Self:
		"""Specify the subset of colours to include."""
		if is_homogenous_type(x, int):
			return super().__new__(cls, (tuple(x),))
		if all(is_homogenous_type(i, int) for i in x):
			return super().__new__(cls, tuple(tuple(i) for i in x))
		errmsg = (
			"SubChromaticInclusion must be initialised with a collection of "
			"integers or collections of integers."
		)
		raise TypeError(errmsg)

	def _simplex_in_domain(
		self,
		column: tuple[list[int], int, float, list[int]],
		filtration: Filtration,  # noqa: ARG002
	) -> bool:
		return any(set(column[3]).issubset(s) for s in self)


class KChromaticInclusion(int, FiltrationInclusion):
	"""
	Corresponds to the inclusion of the :math:`k`-chromatic subcomplex of a chromatic filtration.

	The :math:`k`-chromatic subcomplex is the subfiltration spanned by
	simplices having at most :math:`k` colours.
	"""

	def __new__(cls, k: int) -> Self:
		"""Specify the number of colours :math:`k` for the :math:`k`-chromatic inclusion."""
		if k <= 0:
			errmsg = "KChromaticInclusion must be initialised with a positive integer."
			raise ValueError(errmsg)
		return super().__new__(cls, k)

	def _simplex_in_domain(
		self,
		column: tuple[list[int], int, float, list[int]],
		filtration: Filtration,  # noqa: ARG002
	) -> bool:
		return len(column[3]) <= self


class KChromaticQuotient(int, FiltrationMorphism):
	r"""
	Corresponds to gluing all subfiltrations spanned by :math:`k` colours.

	Given a chromatic filtration :math:`K` with colours :math:`C`,
	this is the quotient map :math:`\bigsqcup_{\substack{I \subset C\\ |I| = k}} K_I \to K`,
	where :math:`K_I` is the subfiltration spanned by simplices
	whose colours are a subset of :math:`I`.
	"""

	def __new__(cls, k: int) -> Self:
		"""Initialise this method."""
		if k <= 0:
			errmsg = "KChromaticQuotient must be initialised with a positive integer."
			raise ValueError(errmsg)
		return super().__new__(cls, k)

	def _compute_diagrams(self, filtration: Filtration, max_dgm_dim: int) -> SixPack:
		BoundaryColumn = NewType("BoundaryColumn", tuple[float, int, list[int]])
		BoundaryMatrix = NewType("BoundaryMatrix", list[BoundaryColumn])

		max_dgm_dim = max(0, filtration.dimension - 1, max_dgm_dim)

		simplex_colours = sorted({vertex.colours[0] for vertex in filtration.simplices[0].values()})
		colour_combinations = sorted(combinations(simplex_colours, self))

		# List of columns in the boundary matrix of the domain and codomain.
		codomain_matrix: BoundaryMatrix = BoundaryMatrix([])
		domain_matrix: BoundaryMatrix = BoundaryMatrix([])
		# List of mappings from domain columns to codomain indices.
		mapping: list[list[int]] = []

		# Offsets[i] is a dictionary that maps a colour subset J
		# to the index of the column in the domain matrix
		# corresponding to the Jth copy of the ith simplex.
		offsets: list[dict[tuple[int, ...], int]] = []

		# Build the matrices
		counter = 0
		for cod_idx, column in enumerate(filtration.serialised()):
			facet_idxs = column[0]
			dimension = max(0, len(facet_idxs) - 1)
			# Only need k-skeleton for k <= max_dgm_dim + 1
			if dimension > max_dgm_dim + 1:
				break
			entrance_time = column[2]
			codomain_matrix.append(BoundaryColumn((entrance_time, dimension, facet_idxs)))
			simplex_colours = column[3]

			# Construct the domain matrix
			offsets.append({})  # Initialize offsets for this codomain column
			for subset in colour_combinations:
				if set(simplex_colours).issubset(subset):
					# If this simplex is J-coloured, we add it to the domain matrix
					# while remembering that for the Jth copy of this simplex i,
					# the index in the domain matrix is counter.
					shifted_facet_idxs = [offsets[facet_idx][subset] for facet_idx in facet_idxs]
					domain_matrix.append(
						BoundaryColumn((entrance_time, dimension, shifted_facet_idxs)),
					)
					mapping.append([cod_idx])
					offsets[cod_idx][subset] = counter
					counter += 1

		logger = logging.getLogger(__name__)
		logger.debug("Domain boundary matrix: \n%s", domain_matrix)
		logger.debug("Codomain boundary matrix: \n%s", codomain_matrix)
		logger.debug("Mappings: \n%s", mapping)

		d, meta = compute_ensemble_cylinder(domain_matrix, codomain_matrix, mapping)

		n_cells_cyl = len(meta.domain) + len(meta.domain_shift) + len(meta.codomain)
		entrance_times = [0.0] * n_cells_cyl
		dimensions = [0] * n_cells_cyl

		for dom_idx, cyl_dom_idx in enumerate(meta.domain):
			entrance_times[cyl_dom_idx] = domain_matrix[dom_idx][0]
			dimensions[cyl_dom_idx] = domain_matrix[dom_idx][1]

		for dom_shift_idx, cyl_dom_shift_idx in enumerate(meta.domain_shift):
			entrance_times[cyl_dom_shift_idx] = domain_matrix[dom_shift_idx][0]
			dimensions[cyl_dom_shift_idx] = domain_matrix[dom_shift_idx][1] + 1

		for cod_idx, cyl_cod_idx in enumerate(meta.codomain):
			entrance_times[cyl_cod_idx] = codomain_matrix[cod_idx][0]
			dimensions[cyl_cod_idx] = codomain_matrix[cod_idx][1]

		return SixPack(
			SimplexPairings(d.ker.paired, d.ker.unpaired),
			SimplexPairings(d.cok.paired, d.cok.unpaired),
			SimplexPairings(d.g.paired, d.g.unpaired),
			SimplexPairings(d.f.paired, d.f.unpaired),
			SimplexPairings(d.im.paired, d.im.unpaired),
			SimplexPairings(d.rel.paired, d.rel.unpaired),
			entrance_times,
			dimensions,
		)


def from_filtration(
	filtration: Filtration,
	mapping_method: FiltrationMorphism,
	max_diagram_dimension: int | None = None,
	threshold: float = 0,
) -> SixPack:
	r"""
	Compute 6-pack of persistence diagrams from a chromatic filtration.

	Given a filtered chromatic simplicial complex :math:`K`
	this function computes the 6-pack of persistence diagram
	associated with a map of :math:`f : L \to K` of filtrations,
	where :math:`L` is some filtration constructed from :math:`K`.

	Args:
		filtration            : A filtered chromatic simplicial complex.
		mapping_method        : The method for constructing the map
			of filtrations.
		max_diagram_dimension :
			Maximum homological dimension for which the persistence diagrams are computed.
			By default diagrams of all dimensions are computed.
		threshold             :
			Retain only points with persistence strictly greater than this value.

	Returns:
		Diagrams corresponding to the following persistence modules
		(where :math:`H_*` is the persistent homology functor and
		:math:`f_*` is the induced map on persistent homology):

		#. :math:`H_*(L)` (domain)
		#. :math:`H_*(K)` (codomain)
		#. :math:`\ker(f_*)` (kernel)
		#. :math:`\mathrm{coker}(f_*)` (cokernel)
		#. :math:`\mathrm{im}(f_*)` (image)
		#. :math:`H_*(K, L)` (relative homology)

		Each diagram is represented by sets of paired and unpaired simplices,
		and contain simplices of all dimensions. ``dgms`` also contains the
		entrance times of the simplices and their dimensions.

	See Also:
		:func:`compute`, :class:`SubChromaticInclusion`,
		:class:`KChromaticInclusion`, :class:`KChromaticQuotient`


	"""
	# If (K, L) is a pair of simplicial complexes where dim(L) <= dim(K) = d
	# then all of the sixpack diagrams for the inclusion map L -> K
	# are zero in dimensions greater than d except for possible the relative diagram.
	# The relative homology is trivial in dimensions greater than d+1.
	if max_diagram_dimension is None:
		max_diagram_dimension = filtration.dimension + 1
	else:
		max_diagram_dimension = min(max_diagram_dimension, filtration.dimension + 1)

	return mapping_method.compute_diagrams(filtration, max_diagram_dimension).threshold(
		max(threshold, 0),
	)


def compute[NumRows: int, NumCols: int](
	points: np.ndarray[tuple[NumRows, NumCols], np.dtype[np.float64]],
	colours: Sequence[int],
	mapping_method: FiltrationMorphism,
	filtration_algorithm: Literal["alpha", "delcech", "delrips"] = "delcech",
	max_diagram_dimension: int | None = None,
	threshold: float = 0,
	max_num_threads: int = 1,
) -> SixPack:
	r"""
	Compute the 6-pack of persistence diagrams of a coloured point-cloud.

	This function constructs a filtered simplicial complex :math:`K`
	from the point cloud, and computes the 6-pack of persistence diagrams
	associated with a map of :math:`f : L \to K` of filtrations,
	where :math:`L` is some filtration constructed from :math:`K`.

	Args:
		points                : Numpy matrix whose columns are points.
		colours               : Sequence of integers describing the colours of the points.
		mapping_method        : The method for constructing the map
			of filtrations.
		filtration_algorithm  :
			Filtration used to construct the chromatic complex.
			Must be one of ``'alpha'``, ``'delcech'``, or ``'delrips'``.
		max_diagram_dimension :
			Maximum homological dimension for which the persistence diagrams are computed.
			By default diagrams of all dimensions are computed.
		threshold             :
			Retain only points with persistence strictly greater than this value.
		max_num_threads       : Maximum number of threads to use to compute the filtration.

	Returns :
		Diagrams corresponding to the following persistence modules
		(where :math:`H_*` is the persistent homology functor and
		:math:`f_*` is the induced map on persistent homology):

		#. :math:`H_*(L)` (domain)
		#. :math:`H_*(K)` (codomain)
		#. :math:`\ker(f_*)` (kernel)
		#. :math:`\mathrm{coker}(f_*)` (cokernel)
		#. :math:`\mathrm{im}(f_*)` (image)
		#. :math:`H_*(K, L)` (relative homology)

		Each diagram is represented by sets of paired and unpaired simplices,
		and contains simplices of all dimensions. ``dgms`` also contains the
		entrance times of the simplices and their dimensions.

	See Also:
		:func:`from_filtration`, :class:`SubChromaticInclusion`,
		:class:`KChromaticInclusion`, :class:`KChromaticQuotient`

	"""
	if isinstance(mapping_method, SubChromaticInclusion) and len(mapping_method) == 1:
		# In this case we are merely checking if the colours of the simplex
		# are in some subset of colours, so we can recolour the points as
		# 0 -> in_subset, 1 -> not_in_subset.
		new_colours: list[int] = (
			np.isin(colours, list(mapping_method), invert=True).astype(int).tolist()
		)
		mapping_method = SubChromaticInclusion((0,))
	else:
		new_colours = list(colours)

	# If X is a d-dimensional point cloud then any subset of X has trivial
	# persistent homology in dimensions greater than or equal to d by a corollary of
	# Alexander duality. This means that in the sixpack of diagrams, only the relative diagram
	# can have non-trivial features in dimension d, and all diagrams are trivial in any dimension
	# greater than d.
	if max_diagram_dimension is None:
		max_diagram_dimension = points.shape[0]
	else:
		max_diagram_dimension = min(max_diagram_dimension, points.shape[0])

	# Compute chromatic complex
	if filtration_algorithm == "delcech":
		filtration, numerical_issues = delcech(points, new_colours, max_num_threads=max_num_threads)
	elif filtration_algorithm == "alpha":
		filtration, numerical_issues = alpha(points, new_colours, max_num_threads=max_num_threads)
	elif filtration_algorithm == "delrips":
		filtration, numerical_issues = delrips(points, new_colours, max_num_threads=max_num_threads)

	if numerical_issues:
		logging.getLogger(__name__).warning(
			"Numerical issues found when computing filtration: %s",
			numerical_issues,
		)
	return from_filtration(
		filtration,
		mapping_method=mapping_method,
		max_diagram_dimension=max_diagram_dimension,
		threshold=threshold,
	)


class SixPack(Mapping):
	"""6-pack of persistence diagrams."""

	__slots__ = ("_dimensions", "_entrance_times", "_simplex_pairings")
	_dimensions: np.ndarray[tuple[int], np.dtype[np.int64]]
	_entrance_times: np.ndarray[tuple[int], np.dtype[np.float64]]
	_simplex_pairings: dict[DiagramName, SimplexPairings]

	type DiagramName = Literal["ker", "cok", "dom", "cod", "im", "rel"]
	"""Names of diagrams in a 6-pack of persistence diagrams."""

	@classmethod
	def names(cls) -> tuple[DiagramName, ...]:
		"""Return the names of the diagrams in the 6-pack."""
		return get_args(cls.DiagramName.__value__)

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

	def threshold(self, tolerance: float) -> SixPack:
		"""Discard all features with persistence ``<=tolerance``."""
		for diagram_name in self:
			pairings = frozenset(
				(b, d)
				for b, d in self[diagram_name].paired
				if self._entrance_times[d] - self._entrance_times[b] > tolerance
			)
			self._simplex_pairings[diagram_name]._paired = pairings  # noqa: SLF001
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
		dimension: list[int] | None = None,
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

		if isinstance(dimension, list):
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
			grp["paired"] = self[dgm_name].paired_as_matrix()

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
		names: list[SixPack.DiagramName] = ["ker", "cok", "dom", "cod", "im", "rel"]
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
		paired: Collection[tuple[int, int]] = set(),
		unpaired: Collection[int] = set(),
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

	def paired_as_matrix(self) -> np.ndarray[tuple[int, Literal[2]], np.dtype[np.int64]]:
		"""Return a matrix representation of the finite persistence features in the diagram."""
		return np.array([[x, y] for (x, y) in self._paired], dtype=np.int64)
