"""Tests for the geometric computations in chalc."""

import itertools
import math

import numpy as np

import chalc as ch


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

	def test_delaunay(self) -> None:
		"""Test chromatic delaunay on a simple example."""
		filtration = ch.chromatic.delaunay(self.points, list(self.colours))
		assert_standard_simplex(3, filtration, is_filtered=False, is_chromatic=True)
		assert filtration.simplices[0][0].colours == 2
		assert all(filtration.simplices[0][i].colours == 1 for i in range(1, 4))
		assert all(filtration.simplices[1][i].colours == 3 for i in range(3))
		assert all(filtration.simplices[1][i].colours == 1 for i in range(3, 6))
		assert all(filtration.simplices[2][i].colours == 3 for i in range(3))
		assert filtration.simplices[3][0].colours == 3

	def test_degenerate_delaunay(self) -> None:
		"""Test the Delaunay construction on degenerate examples."""
		points = np.array([[0, 0], [1, 0], [2, 0]]).T
		colours = [0,0,0]
		delaunay = ch.chromatic.delaunay(points, colours)
		assert delaunay.dimension == 1
		assert len(delaunay.simplices[1]) == 2

		points = np.array([[0, 0], [0, 1], [1, 0], [1, 2]]).T
		colours = [0,0,1,1]
		delaunay = ch.chromatic.delaunay(points, colours)
		assert delaunay.dimension == 2
		assert len(delaunay.simplices[1]) == 5
		assert len(delaunay.simplices[2]) == 2

	def test_alpha(self) -> None:
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

	def test_delcech(self) -> None:
		"""Test chromatic delcech on a simple example."""
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

	def test_delrips(self) -> None:
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

	def test_delcech_is_filtration(self) -> None:
		"""Test that the filtration values are monotonic in the Delcech filtration."""
		rng = np.random.default_rng(self.random_seed)
		dims = [1, 2, 3]
		num_colours = [1, 2, 3]
		for d in dims:
			for s in num_colours:
				points = rng.uniform(size=(d, 200))
				colours = rng.integers(0, s, size=200)
				filtration, numerical_errors = ch.chromatic.delcech(points, list(colours))
				assert not numerical_errors
				assert filtration.is_filtration()

	def test_alpha_is_filtration(self) -> None:
		"""Test that the filtration values are monotonic in the chromatic alpha filtration."""
		rng = np.random.default_rng(self.random_seed)
		dims = [1, 2, 3]
		num_colours = [1, 2, 3]
		for d in dims:
			for s in num_colours:
				points = rng.uniform(size=(d, 200))
				colours = rng.integers(0, s, size=200)
				filtration, numerical_errors = ch.chromatic.alpha(points, list(colours))
				assert not numerical_errors
				assert filtration.is_filtration()

	def test_delrips_is_filtration(self) -> None:
		"""Test that the filtration values are monotonic in the chromatic delrips filtration."""
		rng = np.random.default_rng(self.random_seed)
		dims = [1, 2, 3]
		num_colours = [1, 2, 3]
		for d in dims:
			for s in num_colours:
				points = rng.uniform(size=(d, 200))
				colours = rng.integers(0, s, size=200)
				filtration, numerical_errors = ch.chromatic.delrips(points, list(colours))
				assert not numerical_errors
				assert filtration.is_filtration()

				assert not numerical_errors
				assert filtration.is_filtration()


# class Test_sixpack:
# 	random_seed = 44
# 	def test_k_chromatic(self):
# 		rng = np.random.default_rng(self.random_seed)
# 		points = rng.uniform(size=(2, 1000))
# 		colours = rng.integers(0, 3, size=1000)
# 		dgms = ch.sixpack.compute(points, colours, method='chromatic delcech')


class TestFiltration:
	"""Tests for the filtration module."""

	def test_standard_simplex(self) -> None:
		"""Test the construction of the standard simplex."""
		n = 3
		filtration = ch.filtration.standard_simplex(n)
		assert_standard_simplex(n, filtration, is_filtered=False, is_chromatic=False)

	def test_clique_complex(self) -> None:
		"""Test construction of a clique complex."""
		n = 5
		k = 3
		filtration = ch.filtration.clique_complex(n, k)
		for i in range(k + 1):
			num_simplices = math.comb(n, i + 1)
			assert list(filtration.simplices[i].keys()) == list(range(num_simplices))
			gen_vertices = itertools.combinations(range(n), i + 1)
			for j in range(num_simplices):
				assert filtration.simplices[i][j].dimension == i
				assert filtration.simplices[i][j].vertices == list(next(gen_vertices))
				assert filtration.simplices[i][j].colours == 1
				assert filtration.simplices[i][j].filtration_value == 0
				assert filtration.simplices[i][j].label == j


def assert_standard_simplex(
	n: int,
	filtration: ch.filtration.FilteredComplex,
	*,
	is_filtered: bool,
	is_chromatic: bool,
) -> None:
	"""Assert that the filtration is a standard simplex of dimension n."""
	assert filtration.dimension == n
	assert len(filtration.simplices) == n + 1
	for i in range(n + 1):
		num_simplices = math.comb(n + 1, i + 1)
		assert list(filtration.simplices[i].keys()) == list(range(num_simplices))
		gen_vertices = itertools.combinations(range(n + 1), i + 1)
		for k in range(num_simplices):
			assert filtration.simplices[i][k].dimension == i
			assert filtration.simplices[i][k].vertices == list(next(gen_vertices))
			assert filtration.simplices[i][k].label == k
			if not is_chromatic:
				assert filtration.simplices[i][k].colours == 1
			if not is_filtered:
				assert filtration.simplices[i][k].filtration_value == 0
