"""Test various methods of inclusions in sixpack computation."""

import math
from dataclasses import fields
from itertools import combinations, product
from pathlib import Path

import h5py
import numpy as np
import pytest
from gudhi.bottleneck import bottleneck_distance

import chalc as ch
from chalc.morphisms.types import SixPackMorphisms, SixPackRepresentatives

random_seed = 44


def test_kchromatic_inclusion_vs_quotient() -> None:
	"""Test that k-chromatic inclusion and k-chromatic gluing map yield same result when k=1."""
	rng = np.random.default_rng(random_seed)
	num_points = 500
	dims = [1, 2, 3]
	num_colours = [2, 3]
	for dim, s in product(dims, num_colours):
		points = rng.uniform(size=(dim, num_points))
		colours = rng.integers(0, s, size=num_points).tolist()
		filtration = ch.chromatic.delaunay_rips(points, colours, max_num_threads=0)
		dgms_inclusion = ch.morphisms.KChromaticInclusion(
			filtration,
			1,
		).sixpack()
		dgms_quotient = ch.morphisms.KChromaticQuotient(
			filtration,
			1,
		).sixpack()
		assert dgms_inclusion.max_nonempty_degree() == dgms_quotient.max_nonempty_degree()
		dists = [
			bottleneck_distance(
				getattr(dgms_inclusion, name).intervals(d),
				getattr(dgms_quotient, name).intervals(d),
			)
			for name, _ in dgms_inclusion
			for d in range(dgms_inclusion.max_nonempty_degree() + 1)
		]
		assert math.isclose(max(dists), 0, abs_tol=1e-6)


def test_kchromatic_vs_subchromatic_inclusion() -> None:
	"""Test that k-chromatic and subchromatic inclusions agree when appropriate."""
	rng = np.random.default_rng(random_seed)
	num_points = 500
	dims = [1, 2, 3]
	num_colours = [2, 3]
	for dim, s in product(dims, num_colours):
		for k in range(1, s):
			colour_combos = list(combinations(range(s), k))
			points = rng.uniform(size=(dim, num_points))
			colours = rng.integers(0, s, size=num_points).tolist()
			filtration = ch.chromatic.delaunay_rips(
				points,
				colours,
				max_num_threads=0,
			)
			dgms_kchromatic = ch.morphisms.KChromaticInclusion(
				filtration,
				k,
			).sixpack()
			dgms_subchromatic = ch.morphisms.SubChromaticInclusion(
				filtration,
				colour_combos,
			).sixpack()
			assert dgms_kchromatic.max_nonempty_degree() == dgms_subchromatic.max_nonempty_degree()
			dists = [
				bottleneck_distance(
					getattr(dgms_kchromatic, name).intervals(d),
					getattr(dgms_subchromatic, name).intervals(d),
				)
				for name, _ in dgms_kchromatic
				for d in range(dgms_kchromatic.max_nonempty_degree() + 1)
			]
			assert math.isclose(max(dists), 0, abs_tol=1e-6)


def test_kchromatic_vs_subchromatic_quotient() -> None:
	"""Test that k-chromatic and subchromatic quotients agree when appropriate."""
	rng = np.random.default_rng(random_seed)
	num_points = 500
	dims = [1, 2, 3]
	num_colours = [2, 3]
	for dim, s in product(dims, num_colours):
		for k in range(1, s):
			colour_combos = tuple((c,) for c in combinations(range(s), k))
			points = rng.uniform(size=(dim, num_points))
			colours = rng.integers(0, s, size=num_points).tolist()
			filtration = ch.chromatic.delaunay_rips(
				points,
				colours,
				max_num_threads=0,
			)
			dgms_kchromatic = ch.morphisms.KChromaticQuotient(
				filtration,
				k,
			).sixpack()
			dgms_subchromatic = ch.morphisms.SubChromaticQuotient(
				filtration,
				colour_combos,
			).sixpack()
			assert dgms_kchromatic.max_nonempty_degree() == dgms_subchromatic.max_nonempty_degree()
			dists = [
				bottleneck_distance(
					getattr(dgms_kchromatic, name).intervals(d),
					getattr(dgms_subchromatic, name).intervals(d),
				)
				for name, _ in dgms_kchromatic
				for d in range(dgms_kchromatic.max_nonempty_degree() + 1)
			]
			assert math.isclose(max(dists), 0, abs_tol=1e-6)


def test_subchromatic_inclusion_vs_quotient() -> None:
	"""Test that the subchromatic inclusion and quotient agree for a single subfiltration."""
	rng = np.random.default_rng(random_seed)
	num_points = 500
	dims = [1, 2]
	num_colours = 4
	for dim in dims:
		points = rng.uniform(size=(dim, num_points))
		colours = rng.integers(0, num_colours, size=num_points).tolist()
		# Pick two random faces of Delta^3
		tau = []
		for _i in range(2):
			n_face_verts = rng.integers(1, num_colours)
			face = rng.choice(list(range(num_colours)), size=n_face_verts, replace=False)
			tau.append(face)
		filtration = ch.chromatic.delaunay_rips(
			points,
			colours,
			max_num_threads=0,
		)
		dgms_inclusion = ch.morphisms.SubChromaticInclusion(
			filtration,
			tau,
		).sixpack()
		dgms_quotient = ch.morphisms.SubChromaticQuotient(
			filtration,
			[tau],
		).sixpack()
		assert dgms_inclusion.max_nonempty_degree() == dgms_quotient.max_nonempty_degree()
		dists = [
			bottleneck_distance(
				getattr(dgms_inclusion, name).intervals(d),
				getattr(dgms_quotient, name).intervals(d),
			)
			for name, _ in dgms_inclusion
			for d in range(dgms_inclusion.max_nonempty_degree() + 1)
		]
		assert math.isclose(max(dists), 0, abs_tol=1e-6)


@pytest.mark.parametrize(
	"morphism_type",
	[ch.morphisms.KChromaticInclusion, ch.morphisms.KChromaticQuotient],
)
@pytest.mark.parametrize("get_representatives", [False, True])
def test_sixpack_degree_limit(morphism_type: type, *, get_representatives: bool) -> None:
	filtration = ch.filtration.standard_simplex(2)
	for simplex in filtration:
		simplex.filtration_value = float(simplex.dimension)
	morphism = morphism_type(filtration, 1)
	all_features = morphism.sixpack()
	np.testing.assert_array_equal(all_features.domain.intervals(1), [[1.0, 2.0]])
	limited = morphism.sixpack(0, get_representatives=get_representatives)
	assert limited.max_nonempty_degree() == 0
	for _, module in limited:
		assert all(d <= 0 for d in module.nonempty_degrees())
	assert limited.domain.intervals(1).shape == (0, 2)
	assert all(
		np.all(intervals[:, 0] < intervals[:, 1])
		for _, module in limited
		for intervals in module.intervals().values()
	)
	with pytest.raises(ValueError, match="non-negative"):
		morphism.sixpack(-1)


@pytest.mark.parametrize("get_representatives", [False, True])
def test_quotient_retains_relative_degree_above_filtration(*, get_representatives: bool) -> None:
	filtration = ch.filtration.standard_simplex(0)
	morphism = ch.morphisms.SubChromaticQuotient(filtration, [[{0}], [{0}]])
	pack = morphism.sixpack(get_representatives=get_representatives)
	assert pack.max_nonempty_degree() == 1
	np.testing.assert_array_equal(pack.relative.intervals(1), [[0.0, np.inf]])
	limited = morphism.sixpack(0, get_representatives=get_representatives)
	assert limited.max_nonempty_degree() == 0
	assert limited.relative.intervals(1).shape == (0, 2)
	assert limited.relative.intervals(0).shape == (0, 2)
	if get_representatives:
		assert pack.representatives is not None
		assert limited.representatives is not None
		assert set(pack.representatives.relative) == {0}
		assert limited.representatives.relative == {}
		assert pack.morphisms is not None
		assert pack.morphisms.kernel_to_domain == {0: [0, 1]}


@pytest.mark.parametrize(
	"morphism_type",
	[ch.morphisms.KChromaticInclusion, ch.morphisms.KChromaticQuotient],
)
@pytest.mark.parametrize("max_homological_degree", [None, 0])
def test_sixpack_representatives_and_maps(
	morphism_type: type,
	max_homological_degree: int | None,
	tmp_path: Path,
) -> None:
	morphism = morphism_type(ch.filtration.standard_simplex(2), 1)
	pack = morphism.sixpack(max_homological_degree, get_representatives=True)
	assert isinstance(pack, ch.morphisms.SixPack)
	assert pack == morphism.sixpack(max_homological_degree)
	representatives = pack.representatives
	morphisms = pack.morphisms
	assert isinstance(representatives, SixPackRepresentatives)
	assert isinstance(morphisms, SixPackMorphisms)
	for name, module in pack:
		reps = getattr(representatives, name)
		assert all(module.interval(b) for b in reps)
		assert len(reps) == sum(len(intervals) for intervals in module.intervals().values())
	for field in fields(morphisms):
		source, target = field.name.split("_to_")
		coefficients = getattr(morphisms, field.name)
		assert coefficients.keys() == getattr(representatives, source).keys()
		assert all(
			set(image).issubset(getattr(representatives, target)) for image in coefficients.values()
		)
	with h5py.File(tmp_path / "computed_sixpack.h5", "w") as file:
		pack.to_hdf5(file)
	with h5py.File(tmp_path / "computed_sixpack.h5", "r") as file:
		loaded = ch.morphisms.SixPack.from_hdf5(file)
	assert loaded == pack
	assert loaded.pairings == pack.pairings
	assert loaded.representatives == representatives
	assert loaded.morphisms == morphisms
	assert loaded.metadata is None


def test_kernel_map_retained_at_degree_limit() -> None:
	class VerticesOnlyInclusion(ch.morphisms.FiltrationInclusion):
		def simplex_in_domain(self, column: tuple[list[int], int, float, list[int]]) -> bool:
			return not column[0]

	pack = VerticesOnlyInclusion(ch.filtration.standard_simplex(1)).sixpack(
		0,
		get_representatives=True,
	)
	assert pack.representatives is not None
	assert pack.morphisms is not None
	assert pack.pairings.kernel == {2: None}
	assert pack.kernel.nonempty_degrees() == [0]
	assert pack.representatives.kernel == {0: [0, 1]}
	assert pack.morphisms.kernel_to_domain == {0: [0, 1]}


@pytest.mark.parametrize(
	"morphism_type",
	[ch.morphisms.KChromaticInclusion, ch.morphisms.KChromaticQuotient],
)
def test_sixpack_negative_limit_on_zero_dimensional_filtration(morphism_type: type) -> None:
	morphism = morphism_type(ch.filtration.standard_simplex(0), 1)
	with pytest.raises(ValueError, match="non-negative"):
		morphism.sixpack(-1)
