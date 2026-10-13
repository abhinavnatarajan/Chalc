"""Utilities to compute sheaves of persistence modules from chromatic point clouds."""

from collections.abc import Sequence

import numpy as np

from chalc.chromatic import delaunay_cech
from chalc.filtration import Filtration
from chalc.morphisms import FiltrationInclusion, SubChromaticInclusion

type PersistenceMap = dict[int, list[int]]


class SubfiltrationInclusion(FiltrationInclusion):
	r"""Elementary refinement of a monochromatic Delaunay--Čech filtration.

	Given a bichromatic point cloud :math:`X = X_0 \cup X_1`, this class represents the
	inclusion of the monochromatic Delaunay--Čech filtration
	:math:`\mathrm{Del}\check{C}(X)` into the bichromatic Delaunay--Čech filtration
	:math:`\mathrm{Del}\check{C}(X_0 \cup X_1)`. The inclusion is a homotopy-equivalence
	at each filtration index.
	"""

	__slots__ = ("_subfiltration",)

	def __init__(
		self,
		filtration: Filtration,
		subfiltration: Filtration,
	) -> None:
		"""Specify the point cloud for refinement."""
		self._subfiltration = subfiltration
		super().__init__(filtration)

	def simplex_in_domain(  # ruff: ignore[D102]
		self,
		# column: (facets, label, filtration_time, colours)
		column: tuple[list[int], int, float, list[int]],
	) -> bool:
		facets = column[0]
		dim = len(facets) - 1 if len(facets) > 0 else 0
		label = column[1]
		return label in self._subfiltration.simplices[dim]


class InducedPersistentMorphism:
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
	:math:`X_0 \hookrightarrow X`.
	Therefore, the matrix representations returned by this helper class can be
	used to compose maps of persistence modules across successive inclusions of
	point clouds, using matrix multiplication and zeroing out entries indexed by
	disjoint intervals.

	For computational efficiency, Delaunay--Čech filtrations are used for
	computational efficiency to reduce the number of simplex computations.
	To ensure the compositional consistency of the bases computed by this class, we
	employ the following technique. We first compute the map of persistent homology
	induced by the inclusion of the monochromatic filtration on :math:`X_0` into the
	bichromatic filtration on :math:`X`:

	:math:`i: \mathrm{Del}\check{C}(X_0) \hookrightarrow \mathrm{Del}\check{C}(X_0 \cup X_1)`

	Then we compute the map induced by the inclusion of the monochromatic filtration on
	:math:`X` into the bichromatic filtration on :math:`X`:

	:math:`j: \mathrm{Del}\check{C}(X) \hookrightarrow \mathrm{Del}\check{C}(X_0 \cup X_1)`

	This map is a filtered homotopy-equivalence, and hence induces isomorphisms on
	persistent homology. Finally, we compute the map on persistence induced by
	:math:`j^{-1} \circ i`.
	"""

	def __init__(
		self,
		points: np.ndarray[tuple[int, int], np.dtype[np.floating]],
		colours: Sequence,
	) -> None:
		bichromatic_filt = delaunay_cech(points, colours)
		monochromatic_filt = delaunay_cech(points, [0] * len(colours))
		self.inclusion_i = SubChromaticInclusion(bichromatic_filt, {0})
		self.inclusion_j = SubfiltrationInclusion(bichromatic_filt, monochromatic_filt)

	def induced_morphism(self) -> PersistenceMap:
		morphism_i = self.inclusion_i.sixpack(
			get_representatives=True,
			discard_degenerate_intervals=False,
		)
		morphism_j = self.inclusion_j.sixpack(
			get_representatives=True,
			discard_degenerate_intervals=False,
		)
