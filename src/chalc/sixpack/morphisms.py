"""Helper classes representing morphisms between filtered simplicial complexes."""

from abc import ABC, abstractmethod
from collections.abc import Collection, Iterator, Sized
from itertools import combinations
from typing import Any, TypeIs

import numpy as np
from phimaker import compute_ensemble, compute_ensemble_cylinder

from chalc.filtration import Filtration
from chalc.sixpack.types import SimplexPairings, SixPack

# A face of the simplex of colours.
type ColourSimplexFace = frozenset[int]
# A subcomplex of the simplex of colours, represented by its maximal faces.
type ColourSubcomplex = tuple[ColourSimplexFace, ...]


def _is_homogenous_collection[T](x: Any, typ: type[T]) -> TypeIs[Collection[T]]:
	if isinstance(x, Collection):
		return all(isinstance(i, typ) for i in x)
	return False


class FiltrationMorphism(ABC):
	r"""Abstract map between filtered simplicial complexes.

	The map is a combinatorial proxy for the spatial relationships between
	points in the filtration of different colours.
	This is an abstract class. To specify a morphism, instantiate
	one of its concrete subclasses.

	See Also:
		:class:`SubChromaticInclusion`,
		:class:`KChromaticInclusion`,
		:class:`SubChromaticQuotient`,
		:class:`KChromaticQuotient`.
	"""

	__slots__ = ("filtration",)
	filtration: Filtration

	@abstractmethod
	def sixpack(
		self,
		max_diagram_dimension: int | None = None,
	) -> SixPack:
		r"""
		Compute the 6-pack of persistence diagrams of a coloured point-cloud.

		This function constructs a filtered simplicial complex :math:`K`
		from the point cloud, and computes the 6-pack of persistence diagrams
		associated with a map of :math:`f : L \to K` of filtrations,
		where :math:`L` is some filtration constructed from :math:`K`.

		Args:
			max_diagram_dimension :
				Maximum homological dimension for which the persistence diagrams are computed.
				By default diagrams of all dimensions are computed.

		Returns :
			Diagrams corresponding to the following persistence modules
			(where :math:`H_*` is the persistent homology functor and
			:math:`f_*` is the induced map on persistent homology):

			#. :math:`H_*(L)` (domain)
			#. :math:`H_*(K)` (codomain)
			#. :math:`\ker(f_*)` (kernel)
			#. :math:`\mathrm{coker}(f_*)` (cokernel)
			#. :math:`\mathrm{im}(f_*)` (image)
			#. :math:`H_*(\mathrm{cyl}(f), L)` (relative homology)

			Each diagram is represented by sets of paired and unpaired simplices,
			and contains simplices of all dimensions. ``dgms`` also contains the
			entrance times of the simplices and their dimensions.

		"""

	def __init__(self, filtration: Filtration) -> None:  # noqa: D107
		self.filtration = filtration


class FiltrationInclusion(FiltrationMorphism, ABC):
	"""Abstract inclusion of a filtered subcomplex.

	This class constructs the inclusion of an arbitrary subfiltration
	of a chromatic filtration.

	This is an abstract class. To specify an inclusion map,
	use one of its concrete subclasses or create your own.
	To implement your own inclusion map, you need to implement
	the membership test ``simplex_in_domain`` which checks if a simplex
	is in the domain of the inclusion map.

	See Also:
		:class:`SubChromaticInclusion`, :class:`KChromaticInclusion`.

	"""

	def __init__(self, filtration: Filtration) -> None:  # noqa: D107
		super().__init__(filtration)

	@abstractmethod
	def simplex_in_domain(
		self,
		column: tuple[list[int], int, float, list[int]],
	) -> bool:
		"""Check if a simplex is in the domain of the inclusion map.

		A simplex is identified by its column in the boundary matrix of the filtration.
		This column has the same format as the columns returned by
		:func:`chalc.filtration.Filtration.boundary_matrix`.
		"""

	def sixpack(  # noqa: D102
		self,
		max_diagram_dimension: int | None = None,
	) -> SixPack:
		filtration = self.filtration

		# If f: L \to K is an inclusion map of simplicial complexes where dim(L) <= dim(K) = d
		# then all of the sixpack diagrams for this map are zero in dimensions greater than d
		# except for possibly the relative diagram.
		# The relative homology is trivial in dimensions greater than d+1.
		if max_diagram_dimension is None:
			max_diagram_dimension = filtration.dimension + 1
		else:
			max_diagram_dimension = min(max_diagram_dimension, filtration.dimension + 1)

		codomain_boundary_matrix: list[
			tuple[bool, int, list[int]]
		] = []  # [(in_domain, dimension, facet_idxs)]
		entrance_times: list[float] = []
		dimensions: list[int] = []

		# Mapping from index of a column in the boundary matrix of the filtration
		# to the index of that column in the boundary matrix of the codomain.
		# Not every column in the former ends up in the latter,
		# since we may skip higher dimensional simplices if max_diagram_dimension is set.
		# We initialize with -1 to avoid silent bugs.
		codomain_idx: np.ndarray[tuple[int], np.dtype[np.int64]] = np.array(
			[-1] * len(filtration), dtype=int,
		)

		# Build the codomain boundary matrix
		codomain_idx_counter = 0
		for column_idx, column in enumerate(filtration.boundary_matrix()):
			# Check if we need to skip this column.
			facet_idxs = [codomain_idx[idx] for idx in column[0]]
			dimension = max(0, len(facet_idxs) - 1)
			# H_m(K^n) -> H_m(K) is an isomorphism for m < n.
			# Similarly H_m(L^n) -> H_m(L) is an isomorphism for m < n.
			# Therefore H_m(K^n, L^n) -> H_m(K, L) is an isomorphism by the 5-lemma.
			# If n = max_diagram_dimension + 1, we recover H_m(K), H_m(L), and H_m(K, L)
			# correctly for all m <= max_diagram_dimensions.
			# Consequently we also get the correct kernel, cokernel, and image.
			if dimension > max_diagram_dimension + 1:
				continue
			# If we don't skip this column it gets a codomain index.
			codomain_idx[column_idx] = codomain_idx_counter

			entrance_time = column[2]
			codomain_boundary_matrix.append(
				(
					self.simplex_in_domain(column),
					dimension,
					facet_idxs,
				),
			)
			dimensions.append(dimension)
			entrance_times.append(entrance_time)

			codomain_idx_counter += 1

		d = compute_ensemble(codomain_boundary_matrix)
		return SixPack(
			SimplexPairings(d.ker.paired, d.ker.unpaired),
			SimplexPairings(d.cok.paired, d.cok.unpaired),
			SimplexPairings(d.g.paired, d.g.unpaired),
			SimplexPairings(d.f.paired, d.f.unpaired),
			SimplexPairings(d.im.paired, d.im.unpaired),
			SimplexPairings(d.rel.paired, d.rel.unpaired),
			entrance_times,
			dimensions,
		).threshold(0.0)


class SubChromaticInclusion(FiltrationInclusion, Sized):
	r"""
	Inclusion of a subfiltration spanned by any combination(s) of colours.

	Let :math:`\{0, \ldots, s\}` be a set of colours, let :math:`\Delta^s` be the
	`abstract simplicial complex <https://en.wikipedia.org/wikie/Abstract_simplicial_complex>`_
	whose vertices represent individual colours, and let
	:math:`\tau` be any subcomplex of :math:`\Delta^s`.
	For a filtered simplicial complex :math:`K` on a vertex set :math:`V`,
	and a colouring :math:`\mu:V \to \{0, \ldots, s\}` of :math:`V`,
	let :math:`K/\tau` denote the subfiltration of :math:`K`
	comprising simplices :math:`\sigma \in K` satisfying :math:`\mu(\sigma) \in \tau`.
	This class represents the inclusion map :math:`K/\tau \hookrightarrow K`.

	The complex :math:`\tau` is specified by its maximal faces, or by its maximal face if there
	is only one.

	Examples:
		The inclusion of all monochromatic simplices of colours 0 and 1::

			SubChromaticInclusion(filtration, [[0], [1]]).sixpack()

		The inclusion of any simplex with colours in :math:`\{0, 1\}`,
		(which includes all monochromatic simplices of colours 0 and 1),
		i.e., :math:`\tau = \{\{0, 1\}, \{0\}, \{1\}\}`::

			SubChromaticInclusion(filtration, [[0, 1]]).sixpack()

		In this case since :math:`\tau` has a single maximal face, you can also write the following.
		::

			SubChromaticInclusion(filtration, [0, 1]).sixpack()

		You can also specify more general subsets of colours, for example
		:math:`\tau = \{\{0, 1\}, \{1, 2\}, \{0\}, \{1\}, \{2\}\}`. ::

			SubChromaticInclusion(filtration, [[0, 1], [1, 2]]).sixpack()

	See Also:
		:class:`KChromaticInclusion`, :class:`KChromaticQuotient`, :class:`SubChromaticQuotient`.

	"""

	__slots__ = ("_tau",)
	_tau: ColourSubcomplex

	def __init__(
		self,
		filtration: Filtration,
		tau: Collection[Collection[int]] | Collection[int],
	) -> None:
		"""Specify the subset of colours to include."""
		if _is_homogenous_collection(tau, int):
			# Tau specified as a single maximal face.
			self._tau = (frozenset(tau),)
		elif _is_homogenous_collection(tau, Collection):
			# Tau specified as a collection of maximal faces.
			self._tau = tuple(
				frozenset(  # remove duplicate faces in tau
					frozenset(face) for face in tau
				),
			)
		else:
			errmsg = (
				"SubChromaticInclusion must be initialised with a collection of "
				"integers or collections of integers."
			)
			raise TypeError(errmsg)

		# Check that the colours in tau are valid.
		all_colours = frozenset(
			colour for vertex in filtration.simplices[0].values() for colour in vertex.colours
		)
		if not all(face.issubset(all_colours) for face in self._tau):
			err = "Specified colours in tau are not valid."
			raise ValueError(err)
		super().__init__(filtration)

	def simplex_in_domain(  # noqa: D102
		self,
		column: tuple[list[int], int, float, list[int]],
	) -> bool:
		return any(frozenset(column[3]).issubset(sigma) for sigma in self._tau)

	def __len__(self) -> int:  # noqa: D105
		return len(self._tau)

	def __iter__(self) -> Iterator[frozenset[int]]:  # noqa: D105
		yield from self._tau


class KChromaticInclusion(FiltrationInclusion):
	r"""
	Inclusion of the simplices having at most :math:`k` colours.

	The :math:`k`-chromatic subfiltration is spanned by
	simplices having at most :math:`k` colours. This represents
	a special case of :class:`SubChromaticInclusion`.
	Using the notation from :class:`SubChromaticInclusion`,
	this class corresponds to setting :math:`\tau` to be the
	:math:`k`-skeleton of :math:`\Delta^s`, where
	:math:`\{0, \ldots, s\}` is the set of colours.

	In practical terms, the following code::

		KChromaticInclusion(filtration, k).sixpack()

	should give the same 6-pack of persistence diagrams as this::

		SubChromaticInclusion(
			filtration,
			itertools.combinations(range(n_colours), k)
		).sixpack()

	There is, however, a slight performance benefit to using :class:`KChromaticInclusion`
	over :class:`SubChromaticInclusion` in this situation.

	Examples:
		To consider the inclusion of all monochromatic simplices::

			KChromaticInclusion(filtration, 1).sixpack()

		To consider the inclusion of all simplices spanned by at most two colours::

			KChromaticInclusion(filtration, 2).sixpack()

	See Also:
		:class:`SubChromaticInclusion`, :class:`KChromaticQuotient`, :class:`SubChromaticQuotient`.

	"""

	__slots__ = ("k",)
	k: int
	"""The value k."""

	def __init__(self, filtration: Filtration, k: int) -> None:
		"""Specify the number of colours :math:`k` for the :math:`k`-chromatic inclusion."""
		if not isinstance(k, int) or k <= 0:
			errmsg = "KChromaticInclusion must be initialised with a positive integer."
			raise ValueError(errmsg)
		self.k = k
		super().__init__(filtration)

	def simplex_in_domain(  # noqa: D102
		self,
		column: tuple[list[int], int, float, list[int]],
	) -> bool:
		return len(column[3]) <= self.k


class FiltrationQuotient(FiltrationMorphism, ABC):
	"""Represents the gluing map of a disjoint union of subfiltrations.

	This is an abstract class. To specify a quotient map, use one
	of its concrete subclasses or create your own. To implement
	a quotient map, you need to implement the ``simplex_in_filtration``
	method, and make sure that ``self.num_subfiltrations`` is initialized
	before ``sixpack()`` is called.

	See Also:
		:class:`SubChromaticQuotient`, :class:`KChromaticQuotient`.

	"""

	__slots__ = ("num_subfiltrations",)
	num_subfiltrations: int
	"""The number of subfiltrations in the quotient map."""

	def __init__(self, filtration: Filtration, num_subfiltrations: int) -> None:  # noqa: D107
		self.num_subfiltrations = num_subfiltrations
		super().__init__(filtration)

	@abstractmethod
	def simplex_in_filtration(
		self,
		column: tuple[list[int], int, float, list[int]],
		i: int,
	) -> bool:
		r"""Check if a simplex is in the |ith| subfiltration.

		.. |ith| replace:: i\ :sup:`th`\
		"""

	type BoundaryColumn = tuple[float, int, list[int]]
	type BoundaryMatrix = list[BoundaryColumn]

	def sixpack(  # noqa: D102
		self,
		max_diagram_dimension: int | None = None,
	) -> SixPack:
		filtration = self.filtration

		if max_diagram_dimension is None:
			# If f: L \to K is any map of cell complexes then
			# then the domain, codomain, kernel, cokernel, and image
			# are zero in dimensions greater than max(dim(L), dim(K)).
			# The "relative" homology here is really the relative homology
			# of (cyl(f), L), where cyl(f) is the mapping cylinder of f.
			# This has dimension max(dim(L) + 1, dim(K)).
			# If L is a subcomplex or a disjoint union of subcomplexes of K,
			# then this becomes dim(K) + 1.
			# The relative homology is therefore trivial in dimensions
			# greater than dim(K) + 2.
			max_diagram_dimension = filtration.dimension + 2
		else:
			max_diagram_dimension = min(max_diagram_dimension, filtration.dimension + 2)

		codomain_matrix: FiltrationQuotient.BoundaryMatrix = []
		domain_matrix: FiltrationQuotient.BoundaryMatrix = []
		# List of mappings from domain columns to codomain indices.
		mapping: list[list[int]] = []

		# Mapping from index of a column in the boundary matrix of the filtration
		# to the index of that column in the boundary matrix of the codomain.
		# Not every column in the former ends up in the latter,
		# since we may skip higher dimensional simplices if max_diagram_dimension is set.
		# We initialize with -1 to avoid silent bugs.
		codomain_idx: np.ndarray[tuple[int], np.dtype[np.int64]] = np.array(
			[-1] * len(filtration), dtype=int,
		)

		# For a column with index i in the codomain matrix,
		# offsets[i, j] is its index in the domain matrix
		# corresponding to the copy of the associated simplex
		# belonging to the jth subfiltration.
		# Caution: only some values end up being initialized.
		# We initialize with -1 to avoid silent bugs.
		offsets = np.ones(shape=(len(filtration), self.num_subfiltrations), dtype=int) * -1

		# Build the matrices
		codomain_idx_counter = 0
		domain_idx_counter = 0
		for column_idx, column in enumerate(filtration.boundary_matrix()):
			# Check if we need to skip this column.
			facet_idxs = [codomain_idx[idx] for idx in column[0]]
			dimension = max(0, len(facet_idxs) - 1)
			# H_m(K^n) -> H_m(K) is an iso for m < n.
			# Similarly H_m(L^n) -> H_m(L) is iso for m < n.
			# Therefore H_m(K^n, L^n) -> H_m(K, L) is iso by 5-lemma.
			# If n = max_diagram_dimension + 1, the only need simplices
			# with dimension <= n.
			if dimension > max_diagram_dimension + 1:
				continue

			# If we don't skip this column it gets a codomain index.
			codomain_idx[column_idx] = codomain_idx_counter
			entrance_time = column[2]
			codomain_matrix.append((entrance_time, dimension, facet_idxs))

			# Construct the domain matrix
			for j in range(self.num_subfiltrations):
				if self.simplex_in_filtration(column, j):
					# If this simplex is in the jth subfiltration,
					# we add it to the domain matrix
					# while remembering that for the jth copy of this simplex i,
					# the index in the domain matrix is 'counter'.
					shifted_facet_idxs = [offsets[facet_idx, j] for facet_idx in facet_idxs]
					domain_matrix.append(
						(entrance_time, dimension, shifted_facet_idxs),
					)
					mapping.append([codomain_idx_counter])
					offsets[codomain_idx_counter, j] = domain_idx_counter
					domain_idx_counter += 1

			codomain_idx_counter += 1

		d, meta = compute_ensemble_cylinder(domain_matrix, codomain_matrix, mapping)

		# meta.<domain|codomain|domain_shifted> is a vector
		# mapping from the index of a column in the <domain|codomain>
		# matrix to its index in the mapping cylinder matrix.
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
		).threshold(0.0)


class SubChromaticQuotient(FiltrationQuotient):
	r"""Represents a gluing map in a chromatic filtration.

	Let :math:`\{0, \ldots, s\}` be a set of colours, let :math:`\Delta^s` be the
	`abstract simplicial complex <https://en.wikipedia.org/wikie/Abstract_simplicial_complex>`_
	whose vertices represent individual colours, and let
	:math:`\tau_0, \ldots, \tau_m` be any subcomplexes of :math:`\Delta^s`.
	For a filtered simplicial complex :math:`K` on a vertex set :math:`V`,
	and a colouring :math:`\mu:V \to \{0, \ldots, s\}` of :math:`V`,
	let :math:`K/\tau_i`  denote the subfiltration of :math:`K`
	comprising simplices :math:`\sigma \in K` satisfying :math:`\mu(\sigma) \in \tau_i`
	(:math:`1 \leq i \leq m`). This class represents the quotient map

	.. math::

		\bigsqcup_{i=0}^m K/\tau_i \twoheadrightarrow K,

	Each complex :math:`\tau_i` is specified by its maximal faces, or by its maximal face if there
	is only one.


	Examples:
		If there is only one :math:`\tau_i`, then this is the same as the
		:class:`SubChromaticInclusion` of :math:`\tau_i`.
		For example, both of the following computations produce the same
		6-pack of persistence diagrams, corresponding to
		the inclusion of all monochromatic simplices of colours 0 and 1::

			# Using SubChromaticQuotient
			SubChromaticQuotient(
				filtration,
				[
					[[0, 1]],  # tau_0 = {{0,1}, {0}, {1}}
				],
			).sixpack()

			# Using SubChromaticInclusion
			SubChromaticInclusion(
				filtration,
				[[0,1]],
			).sixpack()

		If the :math:`\tau_i` are disjoint, then this class produces the
		same result as :class:`SubChromaticInclusion`::

			# Using SubChromaticQuotient
			SubChromaticQuotient(
				filtration,
				[
					[
						[0, 1],
					],  # tau_0 = {{0,1}, {0}, {1}}
					[
						[2, 3],
					],  # tau_1 = {{2,3}, {2}, {3}}
				]
			).sixpack()

			# Using SubChromaticInclusion
			SubChromaticInclusion(
				filtration,
				[
					# tau = {{0, 1}, {2, 3}, {0}, {1}, {2}, {3}}
					[0, 1], [2, 3],
				]
			).sixpack()

		In general this is not necessarily the case::

			# Using SubChromaticQuotient - gluing two subfiltrations
			SubChromaticQuotient(
				filtration,
				[
					[
						[0, 1],
					],  # tau_0 = {{0,1}, {0}, {1}}
					[
						[1, 2],
					],  # tau_1 = {{1,2}, {1}, {2}}
				],
			).sixpack()

			# Using SubChromaticInclusion - inclusion of a union of two subfiltrations
			SubChromaticInclusion(
				filtration,
				[
					# tau = {{0, 1}, {1, 2}, {0}, {1}, {2}}
					[0, 1], [1, 2],
				]
			)


	See Also:
		:class:`SubChromaticInclusion`, :class:`KChromaticInclusion`, :class:`KChromaticQuotient`.
	"""

	__slots__ = ("_tau",)
	_tau: tuple[ColourSubcomplex, ...]

	def __init__(
		self,
		filtration: Filtration,
		tau: Collection[Collection[Collection[int]]],
	) -> None:
		"""Initialise this method."""
		self._tau = tuple(
			tuple(frozenset(frozenset(maximal_face) for maximal_face in tau_i)) for tau_i in tau
		)
		# Check if all the colours are valid
		all_colours = frozenset(
			colour for vertex in filtration.simplices[0].values() for colour in vertex.colours
		)
		if not all(face.issubset(all_colours) for tau_i in self._tau for face in tau_i):
			err = "Specified colours in tau are not valid."
			raise ValueError(err)
		num_subfiltrations = len(self._tau)
		super().__init__(filtration, num_subfiltrations)

	def simplex_in_filtration(
		self,
		column: tuple[list[int], int, float, list[int]],
		i: int,
	) -> bool:
		r"""Check if a simplex is in the |ith| subfiltration.

		.. |ith| replace:: i\ :sup:`th`\
		"""
		return any(frozenset(column[3]).issubset(face) for face in self._tau[i])


class KChromaticQuotient(FiltrationQuotient):
	r"""
	Corresponds to gluing all subfiltrations spanned by :math:`k` colours.

	This represents a special case of :class:`SubChromaticQuotient`.
	Using the notation from :class:`SubChromaticQuotient`, this class
	corresponds to having the :math:`\tau_i` range over all combinations of
	:math:`k` colours from the set of colours :math:`\{0, \ldots, s\}`.

	In practical terms, the following code::

		KChromaticInclusion(filtration, k).sixpack()

	should give the same 6-pack of persistence diagrams as this::

		n_colours = len(set(colours))
		SubChromaticQuotient(
			filtration,
			tuple(
				(combination,)
				for combination in combinations(range(n_colours), k))
			)
		)

	Note:
		``KChromaticQuotient(1)`` is essentially the same as
		``KChromaticInclusion(1)`` since both represent the inclusion of all
		monochromatic simplices. You should prefer to use ``KChromaticInclusion(1)``
		for performance reasons, since ``KChromaticQuotient(1)`` will compute
		the mapping cylinder of the inclusion map, which is unnecessary.

	See Also:
		:class:`KChromaticInclusion`, :class:`SubChromaticInclusion`, :class:`SubChromaticQuotient`.

	"""

	__slots__ = ("_tau", "k")
	_tau: tuple[tuple[int, ...], ...]
	k: int
	"""The value k."""

	def __init__(self, filtration: Filtration, k: int) -> None:
		"""Initialise this method."""
		if not isinstance(k, int) or k <= 0:
			errmsg = "KChromaticQuotient must be initialised with a positive integer."
			raise ValueError(errmsg)
		self.k = k
		num_colours = len(
			{colour for vertex in filtration.simplices[0].values() for colour in vertex.colours},
		)
		self._tau = tuple(combinations(range(num_colours), k))
		num_subfiltrations = len(self._tau)
		super().__init__(filtration, num_subfiltrations)

	def simplex_in_filtration(
		self,
		column: tuple[list[int], int, float, list[int]],
		i: int,
	) -> bool:
		r"""Check if a simplex is in the |ith| subfiltration.

		.. |ith| replace:: i\ :sup:`th`\
		"""
		return frozenset(column[3]).issubset(self._tau[i])
