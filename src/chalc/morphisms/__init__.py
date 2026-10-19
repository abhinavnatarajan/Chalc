"""Helper classes representing morphisms between filtered simplicial complexes."""

import logging
from abc import ABC, abstractmethod
from collections.abc import Collection, Iterator, Sequence, Sized
from dataclasses import fields
from itertools import combinations
from typing import TYPE_CHECKING, Any, TypeIs

import numpy as np
import phimaker as pm

from chalc.filtration import Filtration
from chalc.morphisms.types import (
	AugmentedPersistenceModule,
	IndexArray,
	PersistenceModuleName,
	SixPack,
	SixPackMorphisms,
	SixPackPairings,
	SixPackRepresentatives,
)

__all__ = [
	"AugmentedPersistenceModule",
	"FiltrationInclusion",
	"FiltrationMorphism",
	"IndexArray",
	"KChromaticInclusion",
	"KChromaticQuotient",
	"PersistenceModuleName",
	"SixPack",
	"SixPackMorphisms",
	"SixPackPairings",
	"SixPackRepresentatives",
	"SubChromaticInclusion",
	"SubChromaticQuotient",
]

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

	__slots__ = ("_filtration",)
	_filtration: Filtration

	@abstractmethod
	def sixpack(
		self,
		max_homological_degree: int | None = None,
		*,
		get_representatives: bool = False,
	) -> SixPack:
		r"""
		Compute the 6-pack of persistence modules of a coloured point-cloud.

		Compute the modules associated with this map :math:`f : L \to K`,
		where :math:`K` is the supplied filtration and :math:`L` is the domain
		specified by the morphism.

		Args:
			max_homological_degree :
				Maximum homological degree for which the persistence modules are computed.
				By default modules of all dimensions are computed.

		Keyword Args:
			get_representatives:
				If ``True``, also compute representative cycles for each persistence module
				in the 6-pack, and the natural maps between the persistence modules in the 6-pack
				with respect to the interval basis defined by these representative cycles.
				Defaults to ``False``.

		Returns :
			A 6-pack object containing the following persistence modules
			(where :math:`H_*` is the persistent homology functor and
			:math:`f_*` is the induced map on persistent homology):

			#. :math:`H_*(L)` (domain)
			#. :math:`H_*(K)` (codomain)
			#. :math:`\ker(f_*)` (kernel)
			#. :math:`\mathrm{coker}(f_*)` (cokernel)
			#. :math:`\mathrm{im}(f_*)` (image)
			#. :math:`H_*(\mathrm{cyl}(f), L)` (relative homology)

			A :class:`SixPack` exposes its modules as attributes and iterates over
			``(module_name, module)`` pairs. Use :meth:`AugmentedPersistenceModule.intervals`
			on these modules to obtain birth and death times in the desired degrees.
			The original birth and death indices, with ``None`` for infinite deaths,
			are available through :attr:`SixPack.pairings`.

			If ``get_representatives=True``, then the returned
			:class:`SixPack <chalc.morphisms.types.SixPack>` has additional attributes:

			#. ``representatives`` contains representative cycles for each module.
			Each member corresponds to one module, and is a dictionary mapping an
			interval in the module to a cycle representative. Intervals are identified
			by the index of their birth among all interval endpoints in the module.

			#. ``morphisms`` contains the natural maps between the persistence modules
			in the 6-pack. Each persistence map is represented by a binary matrix with
			respect to chosen interval bases of its source and target modules. The interval
			generators are specified by the stored representative cycles, whose classes are
			taken in the quotient defining each module. At each filtration index, the
			representatives of the intervals alive at that index form a basis, and restricting
			the matrix to those source and target intervals gives the linear map at that index.
			For the domain, codomain, image, cokernel, and relative persistence modules,
			the birth time of an interval is the maximum birth time of any simplex in its
			representative cycle. This is not necessarily true for the kernel.
		"""

	def __init__(self, filtration: Filtration) -> None:  # noqa: D107
		self._filtration = filtration


class FiltrationInclusion(FiltrationMorphism, ABC):
	"""Abstract inclusion of a filtered subcomplex.

	This class constructs the inclusion of an arbitrary subfiltration
	of a chromatic filtration.

	This is an abstract class. To specify an inclusion map, use one of its concrete subclasses
	or create your own. To implement your own inclusion map, you need to implement the membership
	test ``simplex_in_domain`` which checks if a simplex is in the domain of the inclusion map.

	Representative cycles, if requested in the call to
	:meth:`sixpack <chalc.morphisms.FiltrationInclusion.sixpack>`, are given as numpy arrays
	of int64 indices; each index identifies a simplex by its codomain boundary columns.
	Representative cycles for the kernel, domain, and image are sums of simplices in the domain,
	while cycle representatives for the codomain and cokernel are sums of codomain simplices.
	The representative cycle for each interval in the relative module is given as a pair of
	chains ``(a, b)`` where ``a=db``.

	See Also:
		:class:`SubChromaticInclusion`, :class:`KChromaticInclusion`.

	"""

	def __init__(self, filtration: Filtration) -> None:  # noqa: D107
		super().__init__(filtration)

	@abstractmethod
	def simplex_in_domain(
		self,
		# column: (facets, label, filtration_time, colours)
		column: tuple[list[int], int, float, list[int]],
	) -> bool:
		"""Check if a simplex is in the domain of the inclusion map.

		A simplex is identified by its column in the boundary matrix of the filtration.
		This column has the same format as the columns returned by
		:func:`chalc.filtration.Filtration.boundary_matrix`.
		"""

	def sixpack(  # noqa: D102
		self,
		max_homological_degree: int | None = None,
		*,
		get_representatives: bool = False,
	) -> SixPack:
		filtration = self._filtration

		if max_homological_degree is not None and max_homological_degree < 0:
			err = (
				f"`max_homological_degree` must be non-negative, but got {max_homological_degree}."
			)
			raise ValueError(err)

		# We do not need all the simplices.
		# Let n = max_homological_degree.
		# H_n(K^{n+2}) -> H_n(K) and H_{n+1}(K^{n+2}) -> H_{n+1}(K) are isomorphisms.
		# Similarly H_n(L^{n+2}) -> H_n(L) and H_{n+1}(L^{n+2}) - H_{n+1}(L) are isomorphisms.
		# Therefore H_{n+1}(K^{n+2}, L^{n+2}) -> H_{n+1}(K, L) is an isomorphism by the 5-lemma
		# applied to the long exact sequence of homology of pairs.
		# Consequently we also get the correct kernel, cokernel, and image.
		# Therefore we only need the boundary matrix of simplices whose
		# dimension is at most max_homological_degree + 2.
		if max_homological_degree is None or max_homological_degree + 2 >= filtration.dimension:
			max_simplex_dimension = -1  # need all simplices
		else:
			max_simplex_dimension = max_homological_degree + 2

		# Feature degrees cannot exceed filtration dimension
		if max_homological_degree is None:
			max_homological_degree = filtration.dimension

		codomain_boundary_matrix: list[list[int]] = []
		entrance_times: list[float] = []
		dimensions: list[int] = []
		cols_in_domain: list[int] = []

		# Build the codomain boundary matrix.
		for idx, column in enumerate(filtration.boundary_matrix(max_simplex_dimension)):
			if self.simplex_in_domain(column):
				cols_in_domain.append(idx)
			facet_idxs = column[0]
			codomain_boundary_matrix.append(facet_idxs)
			dimension = max(0, len(facet_idxs) - 1)
			dimensions.append(dimension)
			entrance_time = column[2]
			entrance_times.append(entrance_time)

		module_ensemble = pm.sixpack_from_inclusion(
			codomain_boundary_matrix,
			dimensions,
			cols_in_domain,
			get_representatives=get_representatives,
			num_threads=int(get_representatives),
			slow=False,
		)

		pairings = SixPackPairings(
			**{
				name: filter_pairings(name, module_ensemble, dimensions, max_homological_degree)
				for name in SixPack.module_names()
			},
		)

		representatives = None
		morphisms = None
		if get_representatives:
			representatives, morphisms = filter_representatives(
				module_ensemble,
				dimensions,
				max_homological_degree,
			)
		return SixPack(
			pairings=pairings,
			entrance_times=entrance_times,
			degrees=dimensions,
			representatives=representatives,
			morphisms=morphisms,
		)


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

		# Check that the colours in tau are valid.
		all_colours = frozenset(vertex.colours[0] for vertex in filtration.simplices[0].values())
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

	should give the same 6-pack of persistence modules as this::

	        SubChromaticInclusion(
	            filtration,
	            itertools.combinations(range(n_colours), k),
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
		# Check that the colours in tau are valid.
		all_colours = frozenset(vertex.colours[0] for vertex in filtration.simplices[0].values())
		if k > len(all_colours):
			errmsg = f"Filtration has {len(all_colours)} colours but k={k} was specified."
			logging.getLogger().warning(errmsg)
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
	before :meth:`sixpack <chalc.morphisms.FiltrationQuotient.sixpack>` is called.

	Persistence computations are reduced to the case of a filtered inclusion
	by constructing a mapping cylinder. Therefore, pairing indices in the
	returned sixpack correspond to mapping cylinder indices, not original indices.

	Representative cycles, if requested in the call to
	:meth:``sixpack <chalc.morphisms.FiltrationInclusion.sixpack>``, are given as numpy arrays
	of int64 indices. These indices identify simplices in either the domain or codomain
	by their index in the corresponding boundary matrix. Representative cycles for the kernel
	and domain are arrays of domain indices, while cycle representatives for the image, codomain
	and cokernel are codomain simplices. The representative cycle for each interval in the relative
	module is given as a pair of chains ``(a, b)`` where ``da= 0`` and ``f(a) = db`` (here ``f`` is
	the quotient map). The indices for ``a`` and ``b`` are domain and codomain indices respectively.

	The sixpack returned by this class contains a ``metadata`` dictionary with the
	following keys:

	#. ``domain_indices`` and ``codomain_indices``: dictionaries that map
	domain and codomain indices to their indices in the mapping cylinder.
	#. ``domain_shift``: dictionary that maps each domain index to the index
	of its degree-shifted copy in the mapping cylinder.
	#. ``domain_index_to_subfiltration_index``: a list that maps each
	index in the boundary matrix of the domain of the quotient map to the
	index of the subfiltration to which this column corresponds.
	#. ``domain_index_to_codomain_index``: a list of integers that completely
	describes the quotient map; the list entry at index ``i`` is the index of
	the codomain simplex that is the image of the ``i``th domain simplex.
	Multiple domain indices map to the same codomain index, one for each
	subfiltration.

	See Also:
		:class:`SubChromaticQuotient`, :class:`KChromaticQuotient`.

	"""

	__slots__ = ("_num_subfiltrations",)
	_num_subfiltrations: int

	def __init__(self, filtration: Filtration, num_subfiltrations: int) -> None:
		self._num_subfiltrations = num_subfiltrations
		super().__init__(filtration)

	@abstractmethod
	def simplex_in_filtration(
		self,
		column: tuple[list[int], int, float, list[int]],
		i: int,
	) -> bool:
		r"""Check if a simplex is in the |ith| subfiltration."""

	type BoundaryMatrix = list[tuple[float, int, list[int]]]

	def sixpack(
		self,
		max_homological_degree: int | None = None,
		*,
		get_representatives: bool = False,
	) -> SixPack:
		filtration = self._filtration
		num_subfiltrations = self._num_subfiltrations

		if max_homological_degree is not None and max_homological_degree < 0:
			err = (
				f"`max_homological_degree` must be non-negative, but got {max_homological_degree}."
			)
			raise ValueError(err)

		# We do not need all the simplices.
		# Let n = max_homological_degree.
		# H_n(K^{n+2}) -> H_n(K) and H_{n+1}(K^{n+2}) -> H_{n+1}(K) are isomorphisms.
		# Similarly H_n(L^{n+2}) -> H_n(L) and H_{n+1}(L^{n+2}) - H_{n+1}(L) are isomorphisms.
		# Therefore H_{n+1}(K^{n+2}, L^{n+2}) -> H_{n+1}(K, L) is an isomorphism by the 5-lemma
		# applied to the long exact sequence of homology of pairs.
		# Consequently we also get the correct kernel, cokernel, and image.
		# Therefore we only need the boundary matrix of simplices whose
		# dimension is at most max_homological_degree + 2.
		if max_homological_degree is None or max_homological_degree + 2 >= filtration.dimension:
			max_simplex_dimension = -1  # need all simplices
		else:
			max_simplex_dimension = max_homological_degree + 2

		# Feature degrees cannot exceed the dimension of the boundary matrix.
		# Here the boundary matrix comes from the mapping cylinder, so it has dimension
		# one higher than the dimension of the filtration.
		if max_homological_degree is None:
			max_homological_degree = filtration.dimension + 1

		codomain_matrix: FiltrationQuotient.BoundaryMatrix = []
		domain_matrix: FiltrationQuotient.BoundaryMatrix = []
		# List of mappings from domain columns to codomain indices.
		mapping: list[list[int]] = []

		# For a column with index i in the codomain matrix,
		# offsets[j][i] is its index in the domain matrix
		# corresponding to the copy of the associated simplex
		# belonging to the jth subfiltration.
		# offsets = np.ones(shape=(len(filtration), num_subfiltrations), dtype=int) * -1
		offsets: dict[int, dict[int, int]] = {i: {} for i in range(num_subfiltrations)}

		# Which subfiltration does the domain column belong to.
		dom_to_subfilt_idx: list[int] = []

		# Build the matrices
		domain_idx_counter = 0
		for codomain_idx, column in enumerate(filtration.boundary_matrix(max_simplex_dimension)):
			facet_idxs = column[0]
			entrance_time = column[2]
			dimension = max(0, len(facet_idxs) - 1)
			codomain_matrix.append((entrance_time, dimension, facet_idxs))

			# Construct the domain matrix
			for i in range(num_subfiltrations):
				if self.simplex_in_filtration(column, i):
					# If this simplex is in the jth subfiltration,
					# we add it to the domain matrix
					# while remembering that for the jth copy of this simplex i,
					# the index in the domain matrix is 'counter'.
					shifted_facet_idxs = [offsets[i][facet_idx] for facet_idx in facet_idxs]
					domain_matrix.append(
						(entrance_time, dimension, shifted_facet_idxs),
					)
					mapping.append([codomain_idx])
					offsets[i][codomain_idx] = domain_idx_counter
					dom_to_subfilt_idx.append(i)
					domain_idx_counter += 1

		module_ensemble, meta = pm.sixpack(
			domain_matrix,
			codomain_matrix,
			mapping,
			get_representatives=get_representatives,
			num_threads=int(get_representatives),
			slow=False,
		)

		pairings = SixPackPairings(
			**{
				name: filter_pairings(name, module_ensemble, meta.degrees, max_homological_degree)
				for name in SixPack.module_names()
			},
		)

		metadata = {
			"domain_indices": meta.domain_indices,
			"codomain_indices": meta.codomain_indices,
			"domain_shift": meta.domain_shift,
			"domain_index_to_subfiltration_index": dom_to_subfilt_idx,
			"domain_index_to_codomain_index": [i[0] for i in mapping],
		}
		representatives = None
		morphisms = None
		if get_representatives:
			representatives, morphisms = filter_representatives(
				module_ensemble,
				meta.degrees,
				max_homological_degree,
			)

		sixpack = SixPack(
			pairings=pairings,
			entrance_times=meta.times,
			degrees=meta.degrees,
			representatives=representatives,
			morphisms=morphisms,
		)
		sixpack.metadata = metadata
		return sixpack


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

	Each complex :math:`\tau_i` is specified as a collection of its maximal faces.


	Examples:
		If there is only one :math:`\tau_i`, then this is the same as the
		:class:`SubChromaticInclusion` of :math:`\tau_i`.
		For example, both of the following computations produce the same
		6-pack of persistence modules, corresponding to
		the inclusion of all simplices spanned by colours 0 and 1::

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
	                    [[0, 1]],
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
	                    ],
	                ).sixpack()

	                # Using SubChromaticInclusion
	                SubChromaticInclusion(
	                    filtration,
	                    [
	                        # tau = {{0, 1}, {2, 3}, {0}, {1}, {2}, {3}}
	                        [0, 1],
	                        [2, 3],
	                    ],
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
	                        [0, 1],
	                        [1, 2],
	                    ],
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
		all_colours = frozenset(vertex.colours[0] for vertex in filtration.simplices[0].values())
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
		r"""Check if a simplex is in the |ith| subfiltration."""
		return any(frozenset(column[3]).issubset(face) for face in self._tau[i])


class KChromaticQuotient(FiltrationQuotient):
	r"""
	Corresponds to gluing all subfiltrations spanned by :math:`k` colours.

	This represents a special case of :class:`SubChromaticQuotient`.
	Using the notation from :class:`SubChromaticQuotient`, this class
	corresponds to having the :math:`\tau_i` range over all combinations of
	:math:`k` colours from the set of colours :math:`\{0, \ldots, s\}`.

	In practical terms, the following code::

	        KChromaticQuotient(filtration, k).sixpack()

	should give the same 6-pack of persistence modules as this::

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

		# Warning if k is too large.
		all_colours = frozenset(vertex.colours[0] for vertex in filtration.simplices[0].values())
		if k > len(all_colours):
			errmsg = f"Filtration has {len(all_colours)} colours but k={k} was specified."
			logging.getLogger().warning(errmsg)

		# Initialise tau with combinations from all colours.
		# We cannot assume that the colours are consecutive or start at zero.
		self._tau = tuple(combinations(all_colours, k))
		num_subfiltrations = len(self._tau)
		super().__init__(filtration, num_subfiltrations)

	def simplex_in_filtration(
		self,
		column: tuple[list[int], int, float, list[int]],
		i: int,
	) -> bool:
		r"""Check if a simplex is in the |ith| subfiltration."""
		return frozenset(column[3]).issubset(self._tau[i])


def filter_representatives(
	ensemble: pm.DiagramEnsemble,
	degrees: Sequence[int],
	max_homological_degree: int,
) -> tuple[SixPackRepresentatives, SixPackMorphisms]:
	"""Restrict representative cycles and natural maps to the requested degrees."""
	representatives_dict: dict[PersistenceModuleName, Any] = {}
	for name in SixPack.module_names():
		if name == "relative":
			continue
		reps = getattr(ensemble, f"{name}_representatives")
		representatives_dict[name] = {
			birth_idx: np.asarray(reps[birth_idx], dtype=np.int64)
			for birth_idx in reps
			if degrees[birth_idx] - (name == "kernel") <= max_homological_degree
		}
	relative_reps = ensemble.relative_representatives
	if TYPE_CHECKING:
		assert relative_reps is not None
	representatives_dict["relative"] = {
		birth_idx: (
			np.asarray(relative_reps[birth_idx][0], dtype=np.int64),
			np.asarray(relative_reps[birth_idx][1], dtype=np.int64),
		)
		for birth_idx in relative_reps
		if degrees[birth_idx] <= max_homological_degree
	}
	representatives = SixPackRepresentatives(**representatives_dict)
	morphisms_dict = {}
	for field in fields(SixPackMorphisms):
		morphism_name = field.name
		morph = getattr(ensemble, morphism_name)
		morphisms_dict[morphism_name] = {
			interval_birth_idx: np.asarray(morph[interval_birth_idx], dtype=np.int64)
			for interval_birth_idx in morph
			if degrees[interval_birth_idx] - (morphism_name.startswith("kernel"))
			<= max_homological_degree
		}
	morphisms = SixPackMorphisms(**morphisms_dict)
	return representatives, morphisms


def filter_pairings(
	name: PersistenceModuleName,
	ensemble: pm.DiagramEnsemble,
	degrees: Sequence[int],
	max_homological_degree: int,
) -> dict[int, int | None]:
	"""Restrict raw pairings, accounting for the kernel's degree shift."""
	return {
		b: d
		for b, d in getattr(ensemble, name).items()
		if degrees[b] - (name == "kernel") <= max_homological_degree
	}
