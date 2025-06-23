"""Tests for the geometric computations in chalc."""

import math

import numpy as np

import chalc as ch

from .utils import assert_standard_simplex


class TestChromatic:
	"""Tests for the chromatic module."""

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

	def test_delaunay_triangulation_deterministic(self) -> None:
		"""Test chromatic delaunay on a simple example."""
		filtration = ch.chromatic.delaunay(self.points, list(self.colours))
		assert_standard_simplex(3, filtration, is_filtered=False, is_chromatic=True)
		assert filtration.simplices[0][0].colours == [1]
		assert all(filtration.simplices[0][i].colours == [0] for i in range(1, 4))
		assert all(filtration.simplices[1][i].colours == [0, 1] for i in range(3))
		assert all(filtration.simplices[1][i].colours == [0] for i in range(3, 6))
		assert all(filtration.simplices[2][i].colours == [0, 1] for i in range(3))
		assert filtration.simplices[3][0].colours == [0, 1]

	def test_delaunay_triangulation_degenerate(self) -> None:
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

	def test_chromatic_alpha_deterministic(self) -> None:
		"""Test chromatic alpha on a simple example."""
		filtration, numerical_errors = ch.chromatic.alpha(self.points, list(self.colours))
		assert not numerical_errors
		assert filtration.is_filtration()
		assert all(filtration.simplices[0][i].filtration_value == 0 for i in range(4))
		assert all(math.isclose(filtration.simplices[1][i].filtration_value, 0.5) for i in range(3))
		assert all(
			math.isclose(filtration.simplices[1][i].filtration_value, np.sqrt(3) / 2)
			for i in range(3, 6)
		)
		assert all(
			math.isclose(filtration.simplices[2][i].filtration_value, np.sqrt(3) / 2)
			for i in range(3)
		)
		assert math.isclose(filtration.simplices[2][3].filtration_value, 1)
		assert math.isclose(filtration.simplices[3][0].filtration_value, 1)

	def test_chromatic_delcech_deterministic(self) -> None:
		"""Test chromatic Delaunay--Cech on a simple example."""
		filtration, numerical_errors = ch.chromatic.delcech(self.points, list(self.colours))
		assert not numerical_errors
		assert filtration.is_filtration()
		assert all(filtration.simplices[0][i].filtration_value == 0 for i in range(4))
		assert all(math.isclose(filtration.simplices[1][i].filtration_value, 0.5) for i in range(3))
		assert all(
			math.isclose(filtration.simplices[1][i].filtration_value, np.sqrt(3) / 2)
			for i in range(3, 6)
		)
		assert all(
			math.isclose(filtration.simplices[2][i].filtration_value, np.sqrt(3) / 2)
			for i in range(3)
		)
		assert math.isclose(filtration.simplices[2][3].filtration_value, 1)
		assert math.isclose(filtration.simplices[3][0].filtration_value, 1)

	def test_chromatic_delrips_deterministic(self) -> None:
		"""Test chromatic delrips on a simple example."""
		filtration, numerical_errors = ch.chromatic.delrips(self.points, list(self.colours))
		assert not numerical_errors
		assert filtration.is_filtration()
		assert all(filtration.simplices[0][i].filtration_value == 0 for i in range(4))
		assert all(math.isclose(filtration.simplices[1][i].filtration_value, 0.5) for i in range(3))
		assert all(
			math.isclose(filtration.simplices[1][i].filtration_value, np.sqrt(3) / 2)
			for i in range(3, 6)
		)
		assert all(
			math.isclose(filtration.simplices[2][i].filtration_value, np.sqrt(3) / 2)
			for i in range(4)
		)
		assert math.isclose(filtration.simplices[3][0].filtration_value, np.sqrt(3) / 2)

	def test_chromatic_delcech_random(self) -> None:
		"""Test the chromatic Delaunay--Cech filtration on random data.

		This tests for monotonicity in the filtration values,
		and for numerical instabilities.
		Both single-threaded and multi-threaded implementations are tested,
		and compared for equality.
		"""
		rng = np.random.default_rng(self.random_seed)
		dims = [1, 2, 3]
		num_colours = [1, 2, 3]
		num_points = 200
		for d in dims:
			for s in num_colours:
				points = rng.uniform(size=(d, num_points))
				colours = rng.integers(0, s, size=num_points, dtype=np.uint64)
				filtration, numerical_errors = ch.chromatic.delcech(
					points,
					colours,
					max_num_threads=1,
				)
				assert not numerical_errors
				assert filtration.is_filtration()
				filtration_mt, numerical_errors_mt = ch.chromatic.delcech(
					points,
					colours,
					max_num_threads=0,
				)
				assert not numerical_errors_mt
				assert filtration_mt.is_filtration()
				assert filtration.serialised() == filtration_mt.serialised()

	def test_chromatic_alpha_random(self) -> None:
		"""Test the chromatic alpha filtration on random data.

		This tests for monotonicity in the filtration values,
		and for numerical instabilities.
		Both single-threaded and multi-threaded implementations are tested,
		and compared for equality.
		"""
		rng = np.random.default_rng(self.random_seed)
		dims = [1, 2, 3]
		num_colours = [1, 2, 3]
		num_points = 200
		for d in dims:
			for s in num_colours:
				points = rng.uniform(size=(d, num_points))
				colours = rng.integers(0, s, size=num_points, dtype=np.uint64)
				filtration_mt, numerical_errors_mt = ch.chromatic.alpha(
					points,
					colours,
					max_num_threads=0,
				)
				assert not numerical_errors_mt
				assert filtration_mt.is_filtration()
				filtration, numerical_errors = ch.chromatic.alpha(
					points,
					colours,
					max_num_threads=1,
				)
				assert not numerical_errors
				assert filtration.is_filtration()
				assert filtration.serialised() == filtration_mt.serialised()

	def test_delrips_random(self) -> None:
		"""Test the chromatic Delaunay--Rips filtration on random data.

		This tests for monotonicity in the filtration values,
		and for numerical instabilities.
		Both single-threaded and multi-threaded implementations are tested,
		and compared for equality.
		"""
		rng = np.random.default_rng(self.random_seed)
		dims = [1, 2, 3]
		num_colours = [1, 2, 3]
		for d in dims:
			for s in num_colours:
				points = rng.uniform(size=(d, 200))
				colours = rng.integers(0, s, size=200, dtype=np.uint64)
				filtration, numerical_errors = ch.chromatic.delrips(
					points,
					colours,
					max_num_threads=1,
				)
				assert not numerical_errors
				assert filtration.is_filtration()
				filtration_mt, numerical_errors_mt = ch.chromatic.delrips(
					points,
					colours,
					max_num_threads=0,
				)
				assert not numerical_errors_mt
				assert filtration_mt.is_filtration()
				assert filtration.serialised() == filtration_mt.serialised()
