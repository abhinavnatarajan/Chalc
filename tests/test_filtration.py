import math
from itertools import combinations

import chalc as ch

from .utils import assert_standard_simplex


def test_standard_simplex() -> None:
	"""Test the construction of the standard simplex."""
	n = 3
	filtration = ch.filtration.standard_simplex(n)
	assert_standard_simplex(n, filtration, is_filtered=False, is_chromatic=False)

def test_complete_complex() -> None:
	"""Test construction of a clique complex."""
	n = 5
	k = 3
	filtration = ch.filtration.complete_complex(n, k)
	for i in range(k + 1):
		num_simplices = math.comb(n, i + 1)
		assert list(filtration.simplices[i].keys()) == list(range(num_simplices))
		gen_vertices = combinations(range(n), i + 1)
		for j in range(num_simplices):
			assert filtration.simplices[i][j].dimension == i
			assert filtration.simplices[i][j].vertices == list(next(gen_vertices))
			assert filtration.simplices[i][j].colours == [0]
			assert filtration.simplices[i][j].filtration_value == 0
			assert filtration.simplices[i][j].label == j
