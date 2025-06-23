"""Routines for computing 6-packs of persistence diagrams."""

from __future__ import annotations

import logging
from collections.abc import (
	Collection,
)
from itertools import combinations
from typing import TYPE_CHECKING, Literal, NewType, TypeVar

import numpy as np
from phimaker import compute_ensemble, compute_ensemble_cylinder

from chalc.chromatic import alpha, delcech, delrips

from ._diagram_ensemble import DiagramEnsemble
from ._simplex_pairings import SimplexPairings

if TYPE_CHECKING:
	from collections.abc import (
		Callable,
		Sequence,
	)

	from chalc.filtration import FilteredComplex

__all__ = [
	"DiagramEnsemble",
	"KChromaticGluingMap",
	"KChromaticInclusion",
	"SimplexPairings",
	"SubsetInclusion",
	"compute",
	"from_filtration",
]

NumRows = TypeVar("NumRows", bound=int)
NumCols = TypeVar("NumCols", bound=int)
Size = TypeVar("Size", bound=tuple[int])


class SubsetInclusion(set[int]):
	"""Corresponds to the inclusion of the subcomplex spanned by a given subset of colours."""

	def __init__(self, x: Collection[int]) -> None:
		"""Initialise this method."""
		if not isinstance(x, Collection) or not all(isinstance(i, int) for i in x):
			errmsg = "SubsetInclusion must be initialised with a collection of integers."
			raise TypeError(errmsg)
		super().__init__(x)


class KChromaticInclusion(int):
	"""
	Corresponds to the inclusion of the :math:`k`-chromatic subcomplex of a chromatic filtration.

	Given a chromatic filtration :math:`K`, its :math:`k`-chromatic subcomplex
	is the subfiltration spanned by simplices having at most :math:`k` colours.

	"""

	def __init__(self, k: int) -> None:
		"""Initialise this method."""
		if k <= 0:
			errmsg = "KChromaticInclusion must be initialised with a positive integer."
			raise ValueError(errmsg)
		super().__init__()


class KChromaticGluingMap(int):
	r"""
	Corresponds to gluing all subfiltrations spanned by :math:`k` colours.

	Given a chromatic filtration :math:`K` with colours :math:`C`,
	this is the gluing map :math:`\bigsqcup_{\substack{I \subset C\\ |I| = k}} K_I \to K`,
	where :math:`K_I` is the subfiltration spanned by simplices
	whose colours are a subset of :math:`I`.

	"""

	def __init__(self, k: int) -> None:
		"""Initialise this method."""
		if k <= 0:
			errmsg = "KChromaticGluingMap must be initialised with a positive integer."
			raise ValueError(errmsg)
		super().__init__()


def _get_diagrams_gluing_map(
	filtration: FilteredComplex,
	k: int,
	max_dgm_dim: int,
	threshold: float = 0,
) -> DiagramEnsemble:
	BoundaryColumn = NewType("BoundaryColumn", tuple[float, int, list[int]])
	BoundaryMatrix = NewType("BoundaryMatrix", list[BoundaryColumn])

	max_dgm_dim = max(0, filtration.dimension - 1, max_dgm_dim)

	simplex_colours = sorted({vertex.colours[0] for vertex in filtration.simplices[0].values()})
	colour_combinations = sorted(combinations(simplex_colours, k))

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
	return dgms.threshold(max(threshold, 0))


def _get_diagrams_inclusion(
	filtration: FilteredComplex,
	check_in_domain: Callable[[Collection[int]], bool],
	max_dgm_dim: int,
	threshold: float = 0,
) -> DiagramEnsemble:
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
		dimensions.append(dimension)
		entrance_times.append(column[2])
		matrix.append((check_in_domain(column[3]), dimension, facet_idxs))

	logger = logging.getLogger(__name__)
	logger.debug("Boundary matrix: \n%s", matrix)

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
	return dgms.threshold(max(threshold, 0))


def from_filtration(
	filtration: FilteredComplex,
	mapping_method: SubsetInclusion | KChromaticInclusion | KChromaticGluingMap,
	max_diagram_dimension: int | None = None,
	threshold: float = 0,
) -> DiagramEnsemble:
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
		:func:`compute`, :class:`SubsetInclusion`,
		:class:`KChromaticInclusion`, :class:`KChromaticGluingMap`


	"""
	# If (K, L) is a pair of simplicial complexes where dim(L) <= dim(K) = d
	# then all of the sixpack diagrams for the inclusion map L -> K
	# are zero in dimensions greater than d except for possible the relative diagram.
	# The relative homology is trivial in dimensions greater than d+1.
	if max_diagram_dimension is None:
		max_diagram_dimension = filtration.dimension + 1
	else:
		max_diagram_dimension = min(max_diagram_dimension, filtration.dimension + 1)

	if isinstance(mapping_method, SubsetInclusion):

		def check_in_domain(c: Collection[int]) -> bool:
			return set(c).issubset(mapping_method)

		return _get_diagrams_inclusion(
			filtration,
			check_in_domain,
			max_dgm_dim=max_diagram_dimension,
			threshold=threshold,
		)

	if isinstance(mapping_method, KChromaticInclusion):

		def check_in_domain(c: Collection[int]) -> bool:
			return len(c) <= mapping_method

		return _get_diagrams_inclusion(
			filtration,
			check_in_domain,
			max_dgm_dim=max_diagram_dimension,
			threshold=threshold,
		)

	if isinstance(mapping_method, KChromaticGluingMap):
		return _get_diagrams_gluing_map(
			filtration,
			k=mapping_method,
			max_dgm_dim=max_diagram_dimension,
			threshold=threshold,
		)

	errmsg = "Invalid mapping method."
	raise TypeError(errmsg)


def compute[NumRows: int, NumCols: int](
	points: np.ndarray[tuple[NumRows, NumCols], np.dtype[np.float64]],
	colours: Sequence[int],
	mapping_method: SubsetInclusion | KChromaticInclusion | KChromaticGluingMap,
	filtration_algorithm: Literal["alpha", "delcech", "delrips"] = "delcech",
	max_diagram_dimension: int | None = None,
	threshold: float = 0,
	max_num_threads: int = 1,
) -> DiagramEnsemble:
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
		:func:`from_filtration`, :class:`SubsetInclusion`,
		:class:`KChromaticInclusion`, :class:`KChromaticGluingMap`

	"""
	if isinstance(mapping_method, SubsetInclusion):
		# If computing using the domain, we only need two colours.
		# new colours: 0 -> domain, 0+1 -> codomain
		new_colours: list[int] = (
			np.isin(colours, list(mapping_method), invert=True).astype(int).tolist()
		)
		mapping_method = SubsetInclusion({0})
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
	else:
		errmsg = "Invalid filtration algorithm."
		raise RuntimeError(errmsg)
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
