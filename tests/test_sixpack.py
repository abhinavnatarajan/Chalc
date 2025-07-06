"""Test various methods of inclusions in sixpack computation."""

import math
from itertools import combinations, product

import numpy as np
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
		filtration = ch.chromatic.delrips(points, colours, max_num_threads=0)
		dgms_inclusion = ch.sixpack.KChromaticInclusion(
			filtration,
			1,
		).sixpack()
		dgms_quotient = ch.sixpack.KChromaticQuotient(
			filtration,
			1,
		).sixpack()
		assert dgms_inclusion.max_nonempty_dimension() == dgms_quotient.max_nonempty_dimension()
		dists = [
			bottleneck_distance(
				dgms_inclusion.get_matrix(name, d),
				dgms_quotient.get_matrix(name, d),
			)
			for name in dgms_inclusion
			for d in range(dgms_inclusion.max_nonempty_dimension() + 1)
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
			filtration = ch.chromatic.delrips(
				points,
				colours,
				max_num_threads=0,
			)
			dgms_kchromatic = ch.sixpack.KChromaticInclusion(
				filtration,
				k,
			).sixpack()
			dgms_subchromatic = ch.sixpack.SubChromaticInclusion(
				filtration,
				colour_combos,
			).sixpack()
			assert (
				dgms_kchromatic.max_nonempty_dimension()
				== dgms_subchromatic.max_nonempty_dimension()
			)
			dists = [
				bottleneck_distance(
					dgms_kchromatic.get_matrix(name, d),
					dgms_subchromatic.get_matrix(name, d),
				)
				for name in dgms_kchromatic
				for d in range(dgms_kchromatic.max_nonempty_dimension() + 1)
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
			filtration = ch.chromatic.delrips(
				points,
				colours,
				max_num_threads=0,
			)
			dgms_kchromatic = ch.sixpack.KChromaticQuotient(
				filtration,
				k,
			).sixpack()
			dgms_subchromatic = ch.sixpack.SubChromaticQuotient(
				filtration,
				colour_combos,
			).sixpack()
			assert (
				dgms_kchromatic.max_nonempty_dimension()
				== dgms_subchromatic.max_nonempty_dimension()
			)
			dists = [
				bottleneck_distance(
					dgms_kchromatic.get_matrix(name, d),
					dgms_subchromatic.get_matrix(name, d),
				)
				for name in dgms_kchromatic
				for d in range(dgms_kchromatic.max_nonempty_dimension() + 1)
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
		filtration = ch.chromatic.delrips(
			points,
			colours,
			max_num_threads=0,
		)
		dgms_inclusion = ch.sixpack.SubChromaticInclusion(
			filtration,
			tau,
		).sixpack()
		dgms_quotient = ch.sixpack.SubChromaticQuotient(
			filtration,
			[tau],
		).sixpack()
		assert dgms_inclusion.max_nonempty_dimension() == dgms_quotient.max_nonempty_dimension()
		dists = [
			bottleneck_distance(
				dgms_inclusion.get_matrix(name, d),
				dgms_quotient.get_matrix(name, d),
			)
			for name in dgms_inclusion
			for d in range(dgms_inclusion.max_nonempty_dimension() + 1)
		]
		assert math.isclose(max(dists), 0, abs_tol=1e-6)
