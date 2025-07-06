"""Tests for the geometric computations in chalc."""

import math

import numpy as np
from gudhi.bottleneck import bottleneck_distance

import chalc as ch

from .utils import assert_clique_complex

# equilateral triangle centred at origin
# plus a point at the origin
points = np.array(
	[
		[0, 0],
		[1, 0],
		[np.cos(2 / 3 * np.pi), np.sin(2 / 3 * np.pi)],
		[np.cos(4 / 3 * np.pi), np.sin(4 / 3 * np.pi)],
	],
).T
colours = (1, 0, 0, 0)
random_seed = 44


def test_delaunay_triangulation_deterministic() -> None:
	"""Test chromatic delaunay on a simple example."""
	filtration = ch.chromatic.delaunay(points, list(colours))
	assert_clique_complex(filtration, n=3, k=3, is_filtered=False, is_chromatic=True)
	assert filtration.simplices[0][0].colours == [1]
	assert all(filtration.simplices[0][i].colours == [0] for i in range(1, 4))
	assert all(sorted(filtration.simplices[1][i].colours) == [0, 1] for i in range(3))
	assert all(sorted(filtration.simplices[1][i].colours) == [0] for i in range(3, 6))
	assert all(sorted(filtration.simplices[2][i].colours) == [0, 1] for i in range(3))
	assert sorted(filtration.simplices[3][0].colours) == [0, 1]


def test_delaunay_triangulation_degenerate() -> None:
	"""Test the Delaunay construction on degenerate examples."""
	points = np.array([[0, 0], [1, 0], [2, 0]]).T
	colours = [0, 0, 0]
	delaunay = ch.chromatic.delaunay(points, colours)
	assert delaunay.dimension == 1
	assert len(delaunay.simplices[1]) == 2

	points = np.array([[0, 0], [0, 1], [1, 0], [1, 2]]).T
	colours = [0, 0, 1, 1]
	delaunay = ch.chromatic.delaunay(points, colours)
	assert delaunay.dimension == 2
	assert len(delaunay.simplices[1]) == 5
	assert len(delaunay.simplices[2]) == 2


def test_chromatic_alpha_deterministic() -> None:
	"""Test chromatic alpha on a simple example."""
	filtration = ch.chromatic.alpha(points, list(colours))
	assert filtration.is_filtration()
	assert all(filtration.simplices[0][i].filtration_value == 0 for i in range(4))
	assert all(math.isclose(filtration.simplices[1][i].filtration_value, 0.5) for i in range(3))
	assert all(
		math.isclose(filtration.simplices[1][i].filtration_value, np.sqrt(3) / 2)
		for i in range(3, 6)
	)
	assert all(
		math.isclose(filtration.simplices[2][i].filtration_value, np.sqrt(3) / 2) for i in range(3)
	)
	assert math.isclose(filtration.simplices[2][3].filtration_value, 1)
	assert math.isclose(filtration.simplices[3][0].filtration_value, 1)


def test_chromatic_delcech_deterministic() -> None:
	"""Test chromatic Delaunay--Cech on a simple example."""
	filtration = ch.chromatic.delcech(points, list(colours))
	assert filtration.is_filtration()
	assert all(filtration.simplices[0][i].filtration_value == 0 for i in range(4))
	assert all(math.isclose(filtration.simplices[1][i].filtration_value, 0.5) for i in range(3))
	assert all(
		math.isclose(filtration.simplices[1][i].filtration_value, np.sqrt(3) / 2)
		for i in range(3, 6)
	)
	assert all(
		math.isclose(filtration.simplices[2][i].filtration_value, np.sqrt(3) / 2) for i in range(3)
	)
	assert math.isclose(filtration.simplices[2][3].filtration_value, 1)
	assert math.isclose(filtration.simplices[3][0].filtration_value, 1)


def test_chromatic_delrips_deterministic() -> None:
	"""Test chromatic delrips on a simple example."""
	filtration = ch.chromatic.delrips(points, list(colours))
	assert filtration.is_filtration()
	assert all(filtration.simplices[0][i].filtration_value == 0 for i in range(4))
	assert all(math.isclose(filtration.simplices[1][i].filtration_value, 0.5) for i in range(3))
	assert all(
		math.isclose(filtration.simplices[1][i].filtration_value, np.sqrt(3) / 2)
		for i in range(3, 6)
	)
	assert all(
		math.isclose(filtration.simplices[2][i].filtration_value, np.sqrt(3) / 2) for i in range(4)
	)
	assert math.isclose(filtration.simplices[3][0].filtration_value, np.sqrt(3) / 2)


def test_chromatic_delcech_random() -> None:
	"""Test the chromatic Delaunay--Cech filtration on random data.

	This tests for monotonicity in the filtration values,
	and for numerical instabilities.
	Both single-threaded and multi-threaded implementations are tested,
	and compared for equality.
	"""
	rng = np.random.default_rng(random_seed)
	dims = [1, 2, 3]
	num_colours = [1, 2, 3]
	num_points = 200
	for d in dims:
		for s in num_colours:
			points = rng.uniform(size=(d, num_points))
			colours = rng.integers(0, s, size=num_points).tolist()
			filtration = ch.chromatic.delcech(
				points,
				colours,
				max_num_threads=1,
			)
			assert filtration.is_filtration()
			filtration_mt = ch.chromatic.delcech(
				points,
				colours,
				max_num_threads=0,
			)
			assert filtration_mt.is_filtration()
			assert filtration.boundary_matrix() == filtration_mt.boundary_matrix()


def test_chromatic_alpha_random() -> None:
	"""Test the chromatic alpha filtration on random data.

	This tests for monotonicity in the filtration values,
	and for numerical instabilities.
	Both single-threaded and multi-threaded implementations are tested,
	and compared for equality.
	"""
	rng = np.random.default_rng(random_seed)
	dims = [1, 2, 3]
	num_colours = [1, 2, 3]
	num_points = 200
	for d in dims:
		for s in num_colours:
			points = rng.uniform(size=(d, num_points))
			colours = rng.integers(0, s, size=num_points).tolist()
			filtration_mt = ch.chromatic.alpha(
				points,
				colours,
				max_num_threads=0,
			)
			assert filtration_mt.is_filtration()
			filtration = ch.chromatic.alpha(
				points,
				colours,
				max_num_threads=1,
			)
			assert filtration.is_filtration()
			assert filtration.boundary_matrix() == filtration_mt.boundary_matrix()


def test_delrips_random() -> None:
	"""Test the chromatic Delaunay--Rips filtration on random data.

	This tests for monotonicity in the filtration values,
	and for numerical instabilities.
	Both single-threaded and multi-threaded implementations are tested,
	and compared for equality.
	"""
	rng = np.random.default_rng(random_seed)
	dims = [1, 2, 3]
	num_colours = [1, 2, 3]
	for d in dims:
		for s in num_colours:
			points = rng.uniform(size=(d, 200))
			colours = rng.integers(0, s, size=200).tolist()
			filtration = ch.chromatic.delrips(
				points,
				colours,
				max_num_threads=1,
			)
			assert filtration.is_filtration()
			filtration_mt  = ch.chromatic.delrips(
				points,
				colours,
				max_num_threads=0,
			)
			assert filtration_mt.is_filtration()
			assert filtration.boundary_matrix() == filtration_mt.boundary_matrix()


def test_alpha_delcech_homotopy_equivalent() -> None:
	"""Test that the 6-pack of diagrams from the delcech and alpha filtrations are identical."""
	rng = np.random.default_rng(random_seed)
	dims = [1, 2]
	num_colours = [2, 3]
	num_points = 500
	for dim in dims:
		for s in num_colours:
			points = rng.uniform(size=(dim, num_points))
			colours = rng.integers(0, s + 1, size=num_points).tolist()
			chromatic_delcech = ch.chromatic.delcech(points, colours)
			chromatic_alpha = ch.chromatic.alpha(points, colours)
			dgms_delcech = ch.sixpack.SubChromaticInclusion(chromatic_delcech, {0}).sixpack()
			dgms_alpha = ch.sixpack.SubChromaticInclusion(chromatic_alpha, {0}).sixpack()
			assert dgms_delcech.max_nonempty_dimension() == dgms_alpha.max_nonempty_dimension()
			dists = [
				bottleneck_distance(
					dgms_delcech.get_matrix("cod", d),
					dgms_alpha.get_matrix("cod", d),
				)
				for d in range(dgms_delcech.max_nonempty_dimension() + 1)
			]
			assert math.isclose(max(dists), 0, abs_tol=1e-6)
