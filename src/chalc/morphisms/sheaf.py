"""Utilities to compute sheaves of persistence modules from chromatic point clouds."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import chain, combinations
from typing import TYPE_CHECKING

import numpy as np

from chalc.chromatic import delaunay_cech
from chalc.morphisms import (
	AugmentedPersistenceModule,
	FiltrationInclusion,
	SixPackMorphisms,
	SubChromaticInclusion,
)
from chalc.morphisms.types import _sparse_rows_equal

if TYPE_CHECKING:
	from collections.abc import Sequence

	from chalc.filtration import Filtration
	from chalc.morphisms.types import IndexArray

type ColourSet = frozenset[int]


@dataclass(slots=True, frozen=True, eq=False)
class PersistenceModuleMorphism:
	"""Sparse map in the interval bases of ``domain`` and ``codomain``.

	``map`` sends module-local source birth indices to one-dimensional int64 arrays of module-local
	target birth indices, with coefficients in F2. Zero columns are omitted.
	At a given time, restrict to source and target intervals alive at that time.
	Nonzero coefficients preserve homological degree and, for source interval
	``[b_s, d_s)`` and target interval ``[b_t, d_t)``, require ``b_t <= b_s < d_t <= d_s``.
	"""

	domain: AugmentedPersistenceModule
	codomain: AugmentedPersistenceModule
	map: dict[int, IndexArray]

	__hash__ = None

	def __eq__(self, other: object) -> bool:
		"""Compare the modules and sparse support arrays by value."""
		if not isinstance(other, PersistenceModuleMorphism):
			return NotImplemented
		return (
			self.domain == other.domain
			and self.codomain == other.codomain
			and _sparse_rows_equal(self.map, other.map)
		)

	def __post_init__(self) -> None:
		"""Validate generator indices, homological degrees, and interval compatibility."""
		if not self.map.keys() <= self.domain.pairings.keys():
			msg = "Map keys should refer to domain intervals."
			raise ValueError(msg)
		if not frozenset(chain.from_iterable(self.map.values())).issubset(
			self.codomain.pairings.keys(),
		):
			msg = "Map values should refer to codomain intervals."
			raise ValueError(msg)
		domain_degrees = self.domain.degrees
		codomain_degrees = self.codomain.degrees
		for source, targets in self.map.items():
			source_birth, source_death = self.domain.interval(source)
			for target_index in targets:
				target = int(target_index)
				if domain_degrees[source] != codomain_degrees[target]:
					msg = f"Map coefficient {source} -> {target} must preserve homological degree."
					raise ValueError(msg)
				target_birth, target_death = self.codomain.interval(target)
				if not target_birth <= source_birth < target_death <= source_death:
					msg = (
						f"Map coefficient {source} -> {target} must satisfy "
						"target birth <= source birth < target death <= source death."
					)
					raise ValueError(msg)

	def __matmul__(self, other: PersistenceModuleMorphism) -> PersistenceModuleMorphism:
		"""Return the composition of composable maps of persistence modules."""
		if self.domain != other.codomain:
			msg = "Incompatible morphisms"
			raise RuntimeError(msg)
		composed = {}
		for src_interval, image in other.map.items():
			composed_image = set()
			for interval in image:
				composed_image ^= set(self.map.get(int(interval), ()))
			image_filtered = tuple(
				target
				for target in composed_image
				if self.codomain.interval(target)[1] > other.domain.interval(src_interval)[0]
			)
			if image_filtered:
				composed[src_interval] = np.fromiter(sorted(image_filtered), dtype=np.int64)
		return PersistenceModuleMorphism(other.domain, self.codomain, composed)

	@classmethod
	def identity(cls, module: AugmentedPersistenceModule) -> PersistenceModuleMorphism:
		"""Identity morphism on a given module."""
		return cls(module, module, {b: np.array([b], dtype=np.int64) for b in module.pairings})

	@classmethod
	def zero(
		cls,
		source: AugmentedPersistenceModule,
		target: AugmentedPersistenceModule,
	) -> PersistenceModuleMorphism:
		"""Return the zero morphism."""
		return cls(source, target, {})


class SubfiltrationInclusion(FiltrationInclusion):
	r"""Inclusion of a subfiltration with the same vertex set as its parent.

	Both filtrations must use the same vertex labels ``0, ..., num_vertices - 1``.
	Filtration times on the subfiltration must agree with those of the parent filtration.
	"""

	__slots__ = ("_subfiltration",)

	def __init__(
		self,
		subfiltration: Filtration,
		filtration: Filtration,
	) -> None:
		"""Specify the subfiltration and its parent.

		Raises:
			ValueError: If the filtrations have different vertex universes.
		"""
		if subfiltration.num_vertices != filtration.num_vertices:
			msg = "Subfiltration and parent filtration must have the same vertex universe."
			raise ValueError(msg)
		self._subfiltration = subfiltration
		super().__init__(filtration)

	def simplex_in_domain(  # noqa: D102
		self,
		# column: (facets, label, filtration_time, colours)
		column: tuple[list[int], int, float, list[int]],
	) -> bool:
		facets = column[0]
		dim = len(facets) - 1 if len(facets) > 0 else 0
		label = column[1]
		simplices = self._subfiltration.simplices
		return dim < len(simplices) and label in simplices[dim]


def _induced_morphism(
	bichromatic_codomain: Filtration,
	monochromatic_codomain: Filtration,
) -> PersistenceModuleMorphism:
	r"""
	Compute the map of persistence modules induced by an inclusion of point clouds.

	If :math:`X = X_0 \cup X_1` is a bi-chromatic point cloud, the inclusion
	:math:`X_0 \hookrightarrow X` induces an inclusion of Čech filtrations
	:math:`\check{C}(X_0) \hookrightarrow \check{C}(X)`, and correspondingly a
	morphism of persistence modules :math:`\mathrm{PH}(X_0) \to \mathrm{PH}(X)`
	over the field :math:`\mathbb{F}_2`. This class is a helper class to compute
	that morphism. It computes a basis of generating intervals for the two
	persistence modules, and the morphism is returned as a binary matrix with
	respect to these bases. The bases computed by this class come with an additional
	guarantee: if :math:`Y` is a superset of :math:`X`, then the basis computed for
	:math:`\mathrm{PH}(X)` for the inclusion :math:`X \hookrightarrow Y` is the
	same as the basis computed for :math:`\mathrm{PH}(X)` for the inclusion
	:math:`X_0 \hookrightarrow X`, provided the order of points in :math:`X`
	is preserved between computations.
	Therefore, the matrix representations returned by this helper class can be
	used to compose maps of persistence modules across successive inclusions of
	point clouds, using matrix multiplication and zeroing out entries indexed by
	disjoint intervals.

	Delaunay--Čech filtrations are used to reduce the number of simplex computations.
	To ensure the compositional consistency of the bases computed by this class, we
	employ the following technique. We first compute the map of persistent homology
	induced by the inclusion of the monochromatic filtration on :math:`X_0` into the
	bichromatic filtration on :math:`X`:

	:math:`i: \mathrm{Del}\check{C}(X_0) \hookrightarrow \mathrm{Del}\check{C}(X_0 \cup X_1)`

	Then we compute the map induced by the inclusion of the monochromatic filtration on
	:math:`X` into the bichromatic filtration on :math:`X`:

	:math:`j: \mathrm{Del}\check{C}(X) \hookrightarrow \mathrm{Del}\check{C}(X_0 \cup X_1)`

	This map is a filtered homotopy-equivalence, and hence induces isomorphisms on
	persistent homology. Finally, we compute the map on persistence :math:`PH(j)^{-1} \circ PH(i)`.
	"""
	inclusion_i = SubChromaticInclusion(bichromatic_codomain, {0})
	inclusion_j = SubfiltrationInclusion(monochromatic_codomain, bichromatic_codomain)
	# Notation:
	# fi : M_1 -> M_2 is the persistent homology morphism induced by
	# DelCech(X_0) -> DelCech(X_0 \cup X_1).
	# fj : M_3 -> M_2 is the persistent homology isomorphism induced by
	# DelCech(X) -> DelCech(X_0 \cup X_1).
	# f = fj^(-1) compose fi
	sixpack_i = inclusion_i.sixpack(
		get_representatives=True,
	)
	sixpack_j = inclusion_j.sixpack(
		get_representatives=True,
	)
	m_1 = sixpack_i.domain
	m_3 = sixpack_j.domain

	if TYPE_CHECKING:
		assert isinstance(sixpack_i.morphisms, SixPackMorphisms)
		assert isinstance(sixpack_j.morphisms, SixPackMorphisms)
	fi = sixpack_i.morphisms.domain_to_codomain
	fj = sixpack_j.morphisms.domain_to_codomain
	f = {}
	for source, image in fi.items():
		birth = m_1.interval(source)[0]
		# Solve in the full interval bases, then discard coefficients
		# between disjoint source and target intervals.
		coefficients = np.fromiter(
			(int(target) for target in _solve(fj, image) if birth < m_3.interval(int(target))[1]),
			dtype=np.int64,
		)
		if coefficients.size:
			f[source] = coefficients
	return PersistenceModuleMorphism(domain=m_1, codomain=m_3, map=f)


def _solve(m: dict[int, IndexArray], b: IndexArray) -> IndexArray:
	r"""Solve a sparse linear equation mx = b over :math:`\mathbb{F}_2`.

	``m`` maps column indices to supports; ``b`` is a support in the row basis.
	Supports may be unsorted. Return the column indices of one solution as a sorted
	int64 array, or raise ``ValueError`` if ``b`` is outside the column span. Empty and
	dependent columns are allowed.
	"""
	basis: dict[int, tuple[set[int], set[int]]] = {}
	for index, support_list in m.items():
		support = {int(row) for row in support_list}
		coefficients = {index}
		while support:
			pivot = max(support)
			if pivot not in basis:
				basis[pivot] = support, coefficients
				break
			basis_column, basis_coefficients = basis[pivot]
			support ^= basis_column
			coefficients ^= basis_coefficients

	residual = {int(row) for row in b}
	solution: set[int] = set()
	while residual:
		pivot = max(residual)
		if pivot not in basis:
			msg = "The right-hand side is outside the column span."
			raise ValueError(msg)
		column, coefficients = basis[pivot]
		residual ^= column
		solution ^= coefficients
	return np.fromiter(sorted(solution), dtype=np.int64)


class ChromaticSheaf:
	"""Sheaf of persistence modules arising from a coloured point cloud."""

	__slots__ = ("_all_colours", "_max_num_colours", "_morphisms", "_stalks")
	# _points: np.ndarray[tuple[int, ...], np.dtype[np.float64]]
	# _colours: np.ndarray[tuple[int], np.dtype[np.int64]]
	# _filtrations: dict[ColourSet, Filtration]
	_stalks: dict[ColourSet, AugmentedPersistenceModule]
	_morphisms: dict[tuple[ColourSet, ColourSet], PersistenceModuleMorphism]
	_all_colours: ColourSet
	_max_num_colours: int

	@property
	def max_num_colours(self) -> int:
		"""Maximum number of colours to consider at a time.

		This determines the height to which we restrict the colour poset.
		"""
		return self._max_num_colours

	@property
	def colours(self) -> ColourSet:
		"""Return all colours in the indexing colour poset."""
		return self._all_colours

	def stalk(self, colours: ColourSet) -> AugmentedPersistenceModule:
		"""Return the stalk persistence module over a set of colours."""
		return self._stalks[colours]

	def morphism(self, source: ColourSet, target: ColourSet) -> PersistenceModuleMorphism:
		"""Return the morphism of persistence modules between two sets of colours.

		The first time a morphism is requested, it has to be constructed by composing
		morphisms across edges in the Hasse diagram. Intermediate compositions are
		cached, so future calls to this function are faster.
		"""
		if not source <= target <= self._all_colours:
			msg = "Need `source <= target <= colours`."
			raise ValueError(msg)
		if not len(target) <= self._max_num_colours:
			msg = f"Target must have at most `max_num_colours={self._max_num_colours}` colours"
			raise ValueError(msg)
		if (source, target) in self._morphisms:
			return self._morphisms[(source, target)]
		# Move from source to target one colour at a time.
		# This isn't necessarily the most efficient method;
		# one could maintain a DAG of computed morphisms and
		# find the shortest path from source to target to minimise
		# the number of compositions. But compositions aren't that
		# expensive.
		colours_to_add = sorted(target - source)
		cur_set = source
		composed = PersistenceModuleMorphism.identity(self._stalks[source])
		self._morphisms[(source, source)] = composed
		for c in colours_to_add:
			new_set = cur_set | {c}
			new_composed = self._morphisms.get((source, new_set))
			if new_composed is not None:
				composed = new_composed
			else:
				composed = self._morphisms[(cur_set, new_set)] @ composed
				self._morphisms[(source, new_set)] = composed
			cur_set = new_set
		return self._morphisms[(source, target)]

	@classmethod
	def from_chromatic_points(
		cls,
		points: np.ndarray[tuple[int, int], np.dtype[np.floating]],
		colours: Sequence[int],
		max_num_colours: int,
	) -> ChromaticSheaf:
		"""Construct a chromatic sheaf from a coloured point cloud."""
		if len(colours) != points.shape[1]:
			msg = "Length of `colours` must equal number of points."
			raise ValueError(msg)
		if max_num_colours < 0:
			msg = "`max_num_colours` must be positive."
			raise ValueError(msg)

		colours_arr = np.array(colours, copy=True)
		all_colours: frozenset[int] = frozenset(colours_arr)
		sheaf = cls()
		sheaf._morphisms = {}
		sheaf._stalks = {}
		sheaf._all_colours = all_colours
		sheaf._max_num_colours = min(max_num_colours, len(all_colours))

		# Initialize the module at the empty set
		sheaf._stalks[frozenset()] = AugmentedPersistenceModule()
		# For each 1 <= n <= max_num_colours:
		# For each set I_n of colours of size n and subset I_{n-1} of I_n,
		# we compute the stalks over I_n, I_{n-1}, and the morphism over
		# I_{n-1} -> I_n.
		for n_colours in range(sheaf._max_num_colours, 0, -1):
			for colour_set in (frozenset(c) for c in combinations(all_colours, n_colours)):
				in_set = np.array([c in colour_set for c in colours_arr])
				x = points[:, in_set]
				monochromatic_codomain = delaunay_cech(
					x,
					[0] * x.shape[1],
				)
				for colour_set_missing_one in (
					frozenset(c) for c in combinations(colour_set, n_colours - 1)
				):
					if not colour_set_missing_one:  # empty set
						module = sheaf._stalks.get(colour_set)
						if module is None:  # for e.g. because max_num_colours = 1
							module = (
								SubChromaticInclusion(monochromatic_codomain, {0})
								.sixpack(get_representatives=True)
								.domain
							)
							sheaf._stalks[colour_set] = module
						zero_module = sheaf._stalks[colour_set_missing_one]
						sheaf._morphisms[(colour_set_missing_one, colour_set)] = (
							PersistenceModuleMorphism.zero(
								zero_module,
								module,
							)
						)
						break
					missing_colour = next(iter(colour_set.difference(colour_set_missing_one)))
					# Assign colour 1 to points with the missing colour, 0 else
					bichromatic_colours = colours_arr[in_set] == missing_colour
					bichromatic_codomain = delaunay_cech(x, bichromatic_colours)
					# TODO(Abhinav): get the monochromatic domain from
					# the bichromatic codomain itself. Might need to change
					# the filtration class in C++ to track the monochromatic
					# simplices during construction.
					m = _induced_morphism(
						bichromatic_codomain,
						monochromatic_codomain,
					)
					if sheaf._stalks.get(colour_set) is None:
						sheaf._stalks[colour_set] = m.codomain
					sheaf._morphisms[(colour_set_missing_one, colour_set)] = m
					sheaf._stalks[colour_set_missing_one] = m.domain

		return sheaf
