import itertools
import math

import chalc as ch
import numpy as np


class Test_chromatic:
	# equilateral triangle centred at origin
	# plus a point at the origin
	points = np.array(
		[
			[0, 0],
			[1, 0],
			[np.cos(2 / 3 * np.pi), np.sin(2 / 3 * np.pi)],
			[np.cos(4 / 3 * np.pi), np.sin(4 / 3 * np.pi)],
		]
	).T
	colours = [1, 0, 0, 0]
	random_seed = 44

	def test_delaunay(self):
		K = ch.chromatic.delaunay(self.points, self.colours)
		assert_standard_simplex(3, K, filtered=False, chromatic=True)
		assert K.simplices[0][0].colours == 2
		assert all(K.simplices[0][i].colours == 1 for i in range(1, 4))
		assert all(K.simplices[1][i].colours == 3 for i in range(3))
		assert all(K.simplices[1][i].colours == 1 for i in range(3, 6))
		assert all(K.simplices[2][i].colours == 3 for i in range(3))
		assert K.simplices[3][0].colours == 3

	def test_alpha(self):
		K, numerical_errors = ch.chromatic.alpha(self.points, self.colours)
		assert not numerical_errors
		assert K.is_filtration()
		assert all(K.simplices[0][i].filtration_value == 0 for i in range(4))
		assert all(math.isclose(K.simplices[1][i].filtration_value, 0.5) for i in range(3))
		assert all(
			math.isclose(K.simplices[1][i].filtration_value, np.sqrt(3) / 2) for i in range(3, 6)
		)
		assert all(
			math.isclose(K.simplices[2][i].filtration_value, np.sqrt(3) / 2) for i in range(3)
		)
		assert math.isclose(K.simplices[2][3].filtration_value, 1)
		assert math.isclose(K.simplices[3][0].filtration_value, 1)

	def test_delcech(self):
		K, numerical_errors = ch.chromatic.delcech(self.points, self.colours)
		assert not numerical_errors
		assert K.is_filtration()
		assert all(K.simplices[0][i].filtration_value == 0 for i in range(4))
		assert all(math.isclose(K.simplices[1][i].filtration_value, 0.5) for i in range(3))
		assert all(
			math.isclose(K.simplices[1][i].filtration_value, np.sqrt(3) / 2) for i in range(3, 6)
		)
		assert all(
			math.isclose(K.simplices[2][i].filtration_value, np.sqrt(3) / 2) for i in range(3)
		)
		assert math.isclose(K.simplices[2][3].filtration_value, 1)
		assert math.isclose(K.simplices[3][0].filtration_value, 1)

	def test_delrips(self):
		K, numerical_errors = ch.chromatic.delrips(self.points, self.colours)
		assert not numerical_errors
		assert K.is_filtration()
		assert all(K.simplices[0][i].filtration_value == 0 for i in range(4))
		assert all(math.isclose(K.simplices[1][i].filtration_value, 0.5) for i in range(3))
		assert all(
			math.isclose(K.simplices[1][i].filtration_value, np.sqrt(3) / 2) for i in range(3, 6)
		)
		assert all(
			math.isclose(K.simplices[2][i].filtration_value, np.sqrt(3) / 2) for i in range(4)
		)
		assert math.isclose(K.simplices[3][0].filtration_value, np.sqrt(3) / 2)

	def test_delcech_is_filtration(self):
		rng = np.random.default_rng(self.random_seed)
		dims = [1, 2, 3]
		num_colours = [1, 2, 3]
		for d in dims:
			for s in num_colours:
				points = rng.uniform(size=(d, 200))
				colours = rng.integers(0, s, size=200)
				K, numerical_errors = ch.chromatic.delcech(points, colours)
				assert not numerical_errors
				assert(K.is_filtration())

	def test_alpha_is_filtration(self):
		rng = np.random.default_rng(self.random_seed)
		dims = [1, 2, 3]
		num_colours = [1, 2, 3]
		for d in dims:
			for s in num_colours:
				points = rng.uniform(size=(d, 200))
				colours = rng.integers(0, s, size=200)
				K, numerical_errors = ch.chromatic.alpha(points, colours)
				assert not numerical_errors
				assert(K.is_filtration())

	def test_delrips_is_filtration(self):
		rng = np.random.default_rng(self.random_seed)
		dims = [1, 2, 3]
		num_colours = [1, 2, 3]
		for d in dims:
			for s in num_colours:
				points = rng.uniform(size=(d, 200))
				colours = rng.integers(0, s, size=200)
				K, numerical_errors = ch.chromatic.delrips(points, colours)
				assert not numerical_errors
				assert(K.is_filtration())

# class Test_sixpack:
# 	random_seed = 44
# 	def test_k_chromatic(self):
# 		rng = np.random.default_rng(self.random_seed)
# 		points = rng.uniform(size=(2, 1000))
# 		colours = rng.integers(0, 3, size=1000)
# 		dgms = ch.sixpack.compute(points, colours, method='chromatic delcech')


class Test_filtration:
	def test_standard_simplex(self):
		n = 3
		K = ch.filtration.standard_simplex(n)
		assert_standard_simplex(n, K, filtered=False, chromatic=False)

	def test_clique_complex(self):
		n = 5
		k = 3
		K = ch.filtration.clique_complex(n, k)
		for i in range(k + 1):
			num_simplices = math.comb(n, i + 1)
			assert list(K.simplices[i].keys()) == list(range(num_simplices))
			gen_vertices = itertools.combinations(range(n), i + 1)
			for j in range(num_simplices):
				assert K.simplices[i][j].dimension == i
				assert K.simplices[i][j].vertices == list(next(gen_vertices))
				assert K.simplices[i][j].colours == 1
				assert K.simplices[i][j].filtration_value == 0
				assert K.simplices[i][j].label == j


def assert_standard_simplex(
	n: int, K: ch.filtration.FilteredComplex, filtered: bool, chromatic: bool
):
	assert K.dimension == n
	assert len(K.simplices) == n + 1
	for i in range(n + 1):
		num_simplices = math.comb(n + 1, i + 1)
		assert list(K.simplices[i].keys()) == list(range(num_simplices))
		gen_vertices = itertools.combinations(range(n + 1), i + 1)
		for k in range(num_simplices):
			assert K.simplices[i][k].dimension == i
			assert K.simplices[i][k].vertices == list(next(gen_vertices))
			assert K.simplices[i][k].label == k
			if not chromatic:
				assert K.simplices[i][k].colours == 1
			if not filtered:
				assert K.simplices[i][k].filtration_value == 0
