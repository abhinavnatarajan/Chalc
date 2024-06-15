from __future__ import annotations

from collections.abc import Collection, Sequence
from dataclasses import dataclass, field, fields
from typing import Annotated, Callable, ClassVar, overload

import numpy as np
from h5py import Dataset, Group
from phimaker import compute_ensemble

from chalc.chromatic import alpha, delcech, delrips
from chalc.filtration import FilteredComplex

from ._utils import interpolate_docstring

__doc__ = "Routines for computing 6-packs of persistence diagrams."
__all__ = [
	"from_filtration",
	"compute",
	"SimplexPairings",
	"DiagramEnsemble",
	"save_diagrams",
	"load_diagrams",
]

ChromaticMethod = {
	"chromatic alpha": alpha,
	"chromatic delcech": delcech,
	"chromatic delrips": delrips,
}

BdMatType = list[tuple[bool, int, list[int]]]


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
		SimplexPairings._fromPhimaker(d.ker),
		SimplexPairings._fromPhimaker(d.cok),
		SimplexPairings._fromPhimaker(d.g),
		SimplexPairings._fromPhimaker(d.f),
		SimplexPairings._fromPhimaker(d.im),
		SimplexPairings._fromPhimaker(d.rel),
		entrance_times,
		dimensions,
	)
	return dgms.threshold(tolerance)


@interpolate_docstring()
def from_filtration(
	K: FilteredComplex,
	dom: Collection[int] | int | None = None,
	k: int | None = None,
	max_diagram_dimension: int = 2,
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
		max_diagram_dimension : Maximum homological dimension for which the persistence diagrams are computed.
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
	if dom is not None and k is None:
		if isinstance(dom, int):
			dom = [
				dom,
			]
		colours_bitmask = _colours_to_bitmask(dom)

		def check_in_domain(b):
			return _colours_are_subset(b, colours_bitmask)
	elif k is not None and dom is None:

		def check_in_domain(b):
			return _num_colours_in_bitmask(b) <= k  # k-chromatic simplex
	elif k is None and dom is None:
		raise RuntimeError("At least one of k or dom must be provided")
	else:
		raise RuntimeError("Only one of k or dom is allowed")
	# if(K, L) is a pair of simplicial complexes where dim(L) <= dim(K) = d
	# then all of the sixpack diagrams except for the relative homology
	# can be non-zero upt dimension d, and relative can be nontrivial upto dim = d+1
	max_diagram_dimension = min(max_diagram_dimension, K.dimension + 1)
	return _get_diagrams(K, check_in_domain, max_diagram_dimension, tolerance)


@interpolate_docstring()
def compute(
	x: Annotated[np.ndarray, np.float64],
	colours: Sequence[int],
	dom: Collection[int] | int | None = None,
	k: int | None = None,
	method: str = "chromatic alpha",
	max_diagram_dimension: int = 2,
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
		method                : Filtration used to construct the chromatic complex. Must be one of ``${str(list(ChromaticMethod.keys()))}``.
		max_diagram_dimension : Maximum homological dimension for which the persistence diagrams are computed.
		tolerance             : Retain only points with persistence strictly greater than this value.

	Returns :
		Diagrams corresponding to the following persistence modules (where :math:`H_*` is the persistent homology functor and :math:`f_*` is the induced map on persistent homology):

		#. :math:`H_*(L)` (domain)
		#. :math:`H_*(K)` (codomain)
		#. :math:`\\ker(f_*)` (kernel)
		#. :math:`\\mathrm{coker}(f_*)` (cokernel)
		#. :math:`\\mathrm{im}(f_*)` (image)
		#. :math:`H_*(K, L)` (relative homology)

		Each diagram is represented by lists of paired and unpaired simplices,
		and contains simplices of all dimensions. ``dgms`` also contains the
		entrance times of the simplices and their dimensions.
	"""
	if dom is not None and k is None:
		# new colours: 1 -> domain, 0 -> codomain
		if isinstance(dom, int):
			dom = [
				dom,
			]
		new_colours = list(np.isin(colours, list(dom)).astype(np.int64))

		def check_in_domain(b: int):
			return b == 2
	elif k is not None and dom is None:

		def check_in_domain(b: int):
			return _num_colours_in_bitmask(b) <= k  # k-chromatic simplex

		new_colours = list(colours)
	elif k is None and dom is None:
		raise RuntimeError("At least one of k or dom must be provided")
	else:
		raise RuntimeError("Only one of k or dom is allowed")
	# can have non-zero homology upto dimension d for everything except rel
	# rel can have non-zero homology in dimension d+1
	max_diagram_dimension = min(max_diagram_dimension, x.shape[0])
	# Compute chromatic complex
	if method in ChromaticMethod.keys():
		K, _ = ChromaticMethod[method](x, new_colours)
	else:
		raise RuntimeError(f"method must be one of {ChromaticMethod.keys()}")
	return _get_diagrams(K, check_in_domain, max_diagram_dimension, tolerance)


@dataclass
class SimplexPairings:
	"""
	Persistence diagram object, represented by a
	list of simplex pairings and a list of unpaired simplices.
	"""

	#: Set of tuples of indices of paired simplices.
	paired: list[tuple[int, int]] = field(default_factory=list)
	#: Set of indices of unpaired simplices.
	unpaired: list[int] = field(default_factory=list)

	@classmethod
	def _fromPhimaker(cls, obj) -> SimplexPairings:
		return cls(list(obj.paired), list(obj.unpaired))

	def __str__(self):
		return f"Paired: {self.paired}\nUnpaired: {self.unpaired}"

	@classmethod
	def _from_matrices(
		cls, p: Annotated[np.ndarray, np.int64], u: Annotated[np.ndarray, np.int64]
	) -> SimplexPairings:
		if p.shape[1] != 2:
			raise ValueError(f"p must be a (m, 2) matrix, but received a matrix of size {p.shape}")
		paired: list[tuple[int, int]] = [(p[i, 0], p[i, 1]) for i in range(p.shape[0])]
		unpaired = list(u)
		return cls(paired, unpaired)

	def paired_as_matrix(self) -> Annotated[np.ndarray, np.int64]:
		return np.concatenate(list(self.paired)).reshape((-1, 2))

	def __eq__(self, other) -> bool:
		if isinstance(other, SimplexPairings):
			return set(self.paired) == set(other.paired) and set(self.unpaired) == set(
				other.unpaired
			)
		else:
			return NotImplemented


@dataclass
class DiagramEnsemble:
	"""
	6-pack of persistence diagrams.
	"""

	# List of matrices
	diagram_names: ClassVar[tuple[str, ...]] = ("ker", "cok", "dom", "cod", "im", "rel")
	#: Kernel
	ker: SimplexPairings = field(default_factory=lambda: SimplexPairings())
	#: Cokernel
	cok: SimplexPairings = field(default_factory=lambda: SimplexPairings())
	#: Domain
	dom: SimplexPairings = field(default_factory=lambda: SimplexPairings())
	#: Codomain
	cod: SimplexPairings = field(default_factory=lambda: SimplexPairings())
	#: Image
	im: SimplexPairings = field(default_factory=lambda: SimplexPairings())
	#: Relative
	rel: SimplexPairings = field(default_factory=lambda: SimplexPairings())
	#: Entrance times of the simplices
	entrance_times: list[float] = field(default_factory=lambda: [])
	#: Dimensions of the simplices
	dimensions: list[int] = field(default_factory=lambda: [])

	def __eq__(self, other) -> bool:
		if isinstance(other, DiagramEnsemble):
			return all(
				getattr(self, f.name) == getattr(other, f.name) for f in fields(DiagramEnsemble)
			)
		else:
			return NotImplemented

	def threshold(self, tolerance: float) -> DiagramEnsemble:
		"""
		Discard all features with persistence ``<=tolerance``.
		"""
		for name in DiagramEnsemble.diagram_names:
			pairings = [
				(b, d)
				for b, d in getattr(self, name).paired
				if self.entrance_times[d] - self.entrance_times[b] > tolerance
			]
			setattr(getattr(self, name), "paired", pairings)
		return self

	@overload
	def get(self, diagram_name: str, dim: int) -> Annotated[np.ndarray, np.float64]: ...

	@overload
	def get(
		self, diagram_name: str, dim: list[int] | None = None
	) -> list[Annotated[np.ndarray, np.float64]]: ...

	@interpolate_docstring({"diagram_names": diagram_names})
	def get(self, diagram_name, dim=None):
		"""
		Get a specific diagram as a matrix of birth and death times.

		Args:
			diagram_name : One of ``${str(diagram_names)}``.
			dim          : Dimension(s) of the diagram desired. If a list is provided then a list of matrices is returned, with the order of matrices respecting the order of entries of `dim`. If `dim` is not provided then the returned matrix will contain persistent features from all homological dimensions from zero to ``max(self.dimensions)``.

		Returns:
			An :math:`m \\times 2` matrix whose rows are a pair of birth and death times, or a list of such matrices.
		"""
		if dim is None:
			dim = list(range(max(set(self.dimensions))))
		if isinstance(dim, int):
			# a p-dimensional homology class is captured by a pairing of (p+1) simplices for kernels
			if diagram_name == "ker":
				dim += 1
			diagram = getattr(self, diagram_name)
			entrance_times = np.array(self.entrance_times)
			paired_list = np.array(
				[pair for pair in diagram.paired if self.dimensions[pair[0]] == dim], dtype=int
			)
			unpaired_list = np.array(
				[sigma for sigma in diagram.unpaired if self.dimensions[sigma] == dim], dtype=int
			)
			paired_matrix = entrance_times[paired_list.flatten()].reshape((-1, 2))
			unpaired_matrix = entrance_times[unpaired_list][:, np.newaxis]
			inf_array = np.array([np.inf] * len(unpaired_list))[:, np.newaxis]
			unpaired_matrix = np.concatenate((unpaired_matrix, inf_array), axis=1)
			return np.concatenate((paired_matrix, unpaired_matrix), axis=0)
		elif isinstance(dim, list):
			return [self.get(diagram_name, d) for d in dim]


def save_diagrams(dgms: DiagramEnsemble, file: Group) -> None:
	"""
	Save a 6-pack of persistence diagrams to a HDF5 file or group.

	Args:
		dgms : 6-pack of diagrams to save to file/group.
		file : A h5py file or group.
	"""
	file["entrance_times"] = np.array(dgms.entrance_times, dtype=np.float64)
	file["dimensions"] = np.array(dgms.dimensions, dtype=np.int64)
	for diagram_name in DiagramEnsemble.diagram_names:
		grp = file.require_group(diagram_name)
		grp["unpaired"] = np.array(list(getattr(dgms, diagram_name).unpaired), dtype=np.int64)
		grp["paired"] = getattr(dgms, diagram_name).paired_as_matrix()


def load_diagrams(file: Group) -> DiagramEnsemble:
	"""
	Load a 6-pack of persistence diagrams from a HDF5 file or group.
	"""
	dgms = DiagramEnsemble()
	if any(f.name not in file for f in fields(DiagramEnsemble)):
		raise RuntimeError("Invalid file: missing fields")

	dgms.entrance_times = (
		list(entrance_times)
		if isinstance((entrance_times := file["entrance_times"]), Dataset)
		and entrance_times.dtype == np.float64
		else []
	)
	dgms.dimensions = (
		list(dimensions)
		if isinstance((dimensions := file["dimensions"]), Dataset) and dimensions.dtype == np.int64
		else []
	)
	for diagram_name in DiagramEnsemble.diagram_names:
		grp = file[diagram_name]
		if isinstance(grp, Group) and all(
			name in list(grp.keys()) for name in ["paired", "unpaired"]
		):
			paired = grp["paired"]
			unpaired = grp["unpaired"]
			if (
				isinstance(paired, Dataset)
				and paired.dtype == np.int64
				and isinstance(unpaired, Dataset)
				and unpaired.dtype == np.int64
			):
				setattr(
					dgms, diagram_name, SimplexPairings._from_matrices(paired[...], unpaired[...])
				)
			else:
				raise RuntimeError("Invalid file: missing fields")
	return dgms
