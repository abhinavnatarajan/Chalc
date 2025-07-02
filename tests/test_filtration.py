import math

import numpy as np
from gudhi.bottleneck import bottleneck_distance

import chalc as ch

from .utils import assert_clique_complex

random_seed = 44


def test_standard_simplex() -> None:
	"""Test the construction of the standard simplex."""
	n = 3
	filtration = ch.filtration.standard_simplex(n)
	assert_clique_complex(filtration, n, n, is_filtered=False, is_chromatic=False)


def test_complete_complex() -> None:
	"""Test construction of a clique complex."""
	n = 5
	k = 3
	filtration = ch.filtration.complete_complex(n, k)
	assert_clique_complex(filtration, n - 1, k, is_filtered=False, is_chromatic=False)


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
			chromatic_delcech, _ = ch.chromatic.delcech(points, colours)
			chromatic_alpha, _ = ch.chromatic.alpha(points, colours)
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
