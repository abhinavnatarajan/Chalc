"""Routines for computing 6-packs of persistence diagrams."""

from __future__ import annotations

from itertools import combinations
from typing import TYPE_CHECKING, Literal, NewType, TypeVar
from warnings import warn

import numpy as np
from phimaker import compute_ensemble, compute_ensemble_cylinder

from chalc.chromatic import alpha, delcech, delrips

from ._diagram_ensemble import DiagramEnsemble
from ._simplex_pairings import SimplexPairings

if TYPE_CHECKING:
	from collections.abc import (
		Callable,
		Collection,
		Sequence,
	)

	from chalc.filtration import FilteredComplex

__all__ = [
	"DiagramEnsemble",
	# "KChromaticGluingMap",
	# "KChromaticInclusion",
	"SimplexPairings",
	# "SubsetInclusion",
	"compute",
	"from_filtration",
]

ChromaticMethod = {
	"alpha": alpha,
	"delcech": delcech,
	"delrips": delrips,
}


NumRows = TypeVar("NumRows", bound=int)
NumCols = TypeVar("NumCols", bound=int)
Size = TypeVar("Size", bound=tuple[int])

# SubsetInclusion = NewType("SubsetInclusion", set[int])
# KChromaticInclusion = NewType("KChromaticInclusion", int)
# KChromaticGluingMap = NewType("KChromaticGluingMap", int)


def _get_diagrams_sheaf(
	filtration: FilteredComplex,
	k: int,
	max_dgm_dim: int,
	threshold: float = 0,
) -> DiagramEnsemble:
	BoundaryColumn = NewType("BoundaryColumn", tuple[float, int, list[int]])
	BoundaryMatrix = NewType("BoundaryMatrix", list[BoundaryColumn])

	all_colours = sorted({vertex.colours[0] for vertex in filtration.simplices[0].values()})
	colour_combinations = sorted(combinations(all_colours, k))

	# List of columns in the boundary matrix of the entire filtration
	codomain_matrix: BoundaryMatrix = BoundaryMatrix([])

	# We need separate copies of the boundary matrix for each colour subset
	# Colour subset J -> list of columns in the boundary matrix indexed by simplices of colours J
	# Each simplex annotated with its corresponding image the in codomain
	domain_matrix_by_subset: dict[tuple[int, ...], list[tuple[BoundaryColumn, list[int]]]] = {
		subset: [] for subset in colour_combinations
	}

	# Build the codomain matrix and the domain matrices for each colour subset
	for idx, column in enumerate(filtration.serialised()):
		facet_idxs = column[0]
		dimension = max(0, len(facet_idxs) - 1)
		# Only need k-skeleton for k <= max_dgm_dim + 1
		if dimension > max_dgm_dim + 1:
			break
		entrance_time = column[2]
		codomain_matrix.append(BoundaryColumn((entrance_time, dimension, facet_idxs)))
		all_colours = column[3]
		subsets = [subset for subset in colour_combinations if set(all_colours).issubset(subset)]
		for subset in subsets:
			domain_matrix_by_subset[subset].append(
				(BoundaryColumn((entrance_time, dimension, facet_idxs)), [idx]),
			)

	# Reindex the duplicated columns in the domain
	shift = 0
	domain_matrix: BoundaryMatrix = BoundaryMatrix([])
	mapping: list[list[int]] = []
	for subset in colour_combinations:
		for col in domain_matrix_by_subset[subset]:
			entrance_time, dimension, original_facet_idxs = col[0]
			facet_idxs = [idx + shift for idx in original_facet_idxs]
			domain_matrix.append(BoundaryColumn((entrance_time, dimension, facet_idxs)))
			mapping.append(col[1])
		shift += len(domain_matrix_by_subset[subset])

	d, meta = compute_ensemble_cylinder(domain_matrix, codomain_matrix, mapping)
	dom_cyl_idx = {j: i for i, j in enumerate(meta.domain)} | {
		j: i for i, j in enumerate(meta.domain_shift)
	}
	cod_cyl_idx = {j: i for i, j in enumerate(meta.codomain)}
	n_cells_cyl = len(meta.domain) + len(meta.domain_shift) + len(meta.codomain)
	entrance_times = [
		domain_matrix[dom_cyl_idx[i]][0] if i in dom_cyl_idx else codomain_matrix[cod_cyl_idx[i]][0]
		for i in range(n_cells_cyl)
	]
	dimensions = [
		domain_matrix[dom_cyl_idx[i]][1] if i in dom_cyl_idx else codomain_matrix[cod_cyl_idx[i]][1]
		for i in range(n_cells_cyl)
	]

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
	dom: Collection[int] | int | None = None,
	k: int | None = None,
	max_diagram_dimension: int | None = None,
	threshold: float = 0,
) -> DiagramEnsemble:
	r"""
	Compute 6-pack of persistence diagrams from a chromatic filtration.

	Given a filtered chromatic simplicial complex :math:`K`
	and a subcomplex :math:`L` of :math:`K`,
	this function computes the 6-pack of persistence diagram
	associated with the inclusion map :math:`f : L \hookrightarrow K`.
	The subcomplex is specified by the colours of its vertices,
	or by an integer :math:`k` wherein all simplices with
	:math:`k` or fewer colours are considered part of the subcomplex.

	Args:
		filtration            : A filtered chromatic simplicial complex.
		dom                   :
			Integer or collection of integers describing the
			colours of the points in the domain (the subcomplex :math:`L`).
		k                     :
			If not ``None``, then the domain is taken to be the
			:math:`k`-chromatic subcomplex of :math:`K`, i.e.,
			the subcomplex of simplices having at most :math:`k` colours.
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

	"""
	if dom is not None and k is None:
		if isinstance(dom, int):
			dom = [
				dom,
			]

		def check_in_domain(c: Collection[int]) -> bool:
			return set(c).issubset(dom)
	elif k is not None and dom is None:

		def check_in_domain(c: Collection[int]) -> bool:
			return len(c) <= k
	elif k is None and dom is None:
		errmsg = "At least one of k or dom must be provided"
		raise RuntimeError(errmsg)
	else:
		errmsg = "Only one of k or dom is allowed"
		raise RuntimeError(errmsg)

	# If (K, L) is a pair of simplicial complexes where dim(L) <= dim(K) = d
	# then all of the sixpack diagrams for the inclusion map L -> K
	# are zero in dimensions greater than d except for possible the relative diagram.
	# The relative homology is trivial in dimensions greater than d+1.
	if max_diagram_dimension is None:
		max_diagram_dimension = filtration.dimension + 1
	else:
		max_diagram_dimension = min(max_diagram_dimension, filtration.dimension + 1)
	return _get_diagrams_inclusion(filtration, check_in_domain, max_diagram_dimension, threshold)


def compute[NumRows: int, NumCols: int](
	x: np.ndarray[tuple[NumRows, NumCols], np.dtype[np.floating]],
	colours: Sequence[int],
	dom: Collection[int] | int | None = None,
	k: int | None = None,
	method: Literal["alpha", "delcech", "delrips"] = "alpha",
	max_diagram_dimension: int | None = None,
	threshold: float = 0,
) -> DiagramEnsemble:
	r"""
	Compute the 6-pack of persistence diagrams of a coloured point-cloud.

	This function constructs a filtered simplicial complex :math:`K`
	from the point cloud, and computes the 6-pack of persistence diagrams
	associated with the inclusion :math:`f : L \hookrightarrow K`
	where :math:`L` is some filtered subcomplex of :math:`K`.

	Args:
		x                     : Numpy matrix whose columns are points.
		colours               : Sequence of integers describing the colours of the points.
		dom                   :
			Integer or collection of integers describing the colours of
			the points in the domain (the subcomplex :math:`L`).
		k                     :
			If not ``None``, then the domain is taken to be the
			:math:`k`-chromatic subcomplex of :math:`K`, i.e.,
			the subcomplex of simplices having at most :math:`k` colours.
		method                :
			Filtration used to construct the chromatic complex.
			Must be one of ``'alpha'``, ``'delcech'``, or ``'delrips'``.
		max_diagram_dimension :
			Maximum homological dimension for which the persistence diagrams are computed.
			By default diagrams of all dimensions are computed.
		threshold             :
			Retain only points with persistence strictly greater than this value.

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

	"""
	if dom is not None and k is None:
		# If computing using the domain, we only need two colours.
		# new colours: 0 -> domain, 0+1 -> codomain
		if isinstance(dom, int):
			dom = [
				dom,
			]
		new_colours = list(np.isin(colours, list(dom), invert=True).astype(int))
		dom = [0]
	else:
		new_colours = colours

	# If X is a d-dimensional point cloud then any subset of X has trivial
	# persistent homology in dimensions greater than or equal to d by a corollary of
	# Alexander duality. This means that in the sixpack of diagrams, only the relative diagram
	# can have non-trivial features in dimension d, and all diagrams are trivial in any dimension
	# greater than d.
	if max_diagram_dimension is None:
		max_diagram_dimension = x.shape[0]
	else:
		max_diagram_dimension = min(max_diagram_dimension, x.shape[0])

	# Compute chromatic complex
	if method in ChromaticMethod:
		filtration, numerical_issues = ChromaticMethod[method](x, new_colours)
		if numerical_issues:
			warn(
				"Numerical issues found when computing filtration",
				category=RuntimeWarning,
				stacklevel=1,  # warn about
			)
	else:
		errmsg = f"method must be one of {ChromaticMethod.keys()}"
		raise RuntimeError(errmsg)
	return from_filtration(
		filtration,
		dom=dom,
		k=k,
		max_diagram_dimension=max_diagram_dimension,
		threshold=threshold,
	)
