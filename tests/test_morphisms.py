"""Test various methods of inclusions in sixpack computation."""

import math
from itertools import combinations, product

import numpy as np
import pytest
from gudhi.bottleneck import bottleneck_distance

import chalc as ch

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
				dgms_inclusion.get_module(name).intervals(d),
				dgms_quotient.get_module(name).intervals(d),
			)
			for name in dgms_inclusion
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
					dgms_kchromatic.get_module(name).intervals(d),
					dgms_subchromatic.get_module(name).intervals(d),
				)
				for name in dgms_kchromatic
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
					dgms_kchromatic.get_module(name).intervals(d),
					dgms_subchromatic.get_module(name).intervals(d),
				)
				for name in dgms_kchromatic
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
				dgms_inclusion.get_module(name).intervals(d),
				dgms_quotient.get_module(name).intervals(d),
			)
			for name in dgms_inclusion
			for d in range(dgms_inclusion.max_nonempty_degree() + 1)
		]
		assert math.isclose(max(dists), 0, abs_tol=1e-6)


@pytest.mark.parametrize(
	"morphism_type",
	[ch.morphisms.KChromaticInclusion, ch.morphisms.KChromaticQuotient],
)
@pytest.mark.parametrize("discard_degenerate", [False, True])
def test_sixpack_degree_limit(morphism_type: type, *, discard_degenerate: bool) -> None:
	filtration = ch.filtration.standard_simplex(2)
	morphism = morphism_type(filtration, 1)
	all_features = morphism.sixpack(discard_degenerate_intervals=False)
	assert all_features.max_nonempty_degree() >= 1
	limited = morphism.sixpack(0, discard_degenerate_intervals=discard_degenerate)
	assert limited.max_nonempty_degree() == 0
	for name in limited:
		module = limited.get_module(name)
		assert all(d <= 0 for d in module.degrees)
		for intervals in module.intervals().values():
			if discard_degenerate:
				assert np.all(intervals[:, 1] > intervals[:, 0])
	if not discard_degenerate:
		assert any(
			np.any(mat[:, 0] == mat[:, 1])
			for mat in limited.get_module("domain").intervals().values()
		)
	with pytest.raises(ValueError, match="non-negative"):
		morphism.sixpack(-1)


def test_quotient_retains_relative_degree_above_filtration() -> None:
	filtration = ch.filtration.standard_simplex(0)
	morphism = ch.morphisms.SubChromaticQuotient(filtration, [[{0}], [{0}]])
	pack = morphism.sixpack()
	assert pack.max_nonempty_degree() == 1
	np.testing.assert_array_equal(pack.get_module("relative").intervals(1), [[0.0, np.inf]])
	assert not morphism.sixpack(0).get_module("relative")


@pytest.mark.parametrize(
	"morphism_type",
	[ch.morphisms.KChromaticInclusion, ch.morphisms.KChromaticQuotient],
)
def test_sixpack_negative_limit_on_zero_dimensional_filtration(morphism_type: type) -> None:
	morphism = morphism_type(ch.filtration.standard_simplex(0), 1)
	with pytest.raises(ValueError, match="non-negative"):
		morphism.sixpack(-1)
