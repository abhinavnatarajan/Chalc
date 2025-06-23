"""Test various methods of inclusions in sixpack computation."""

import math
from itertools import product
from typing import get_args

import numpy as np
from gudhi.bottleneck import bottleneck_distance

import chalc as ch


class TestSixpack:
	"""Test the computation of sixpacks."""

	random_seed = 44

	def test_kchromatic_vs_gluing(self) -> None:
		"""Test that k-chromatic inclusion and k-chromatic gluing map yield same result when k=1."""
		rng = np.random.default_rng(self.random_seed)
		num_points = 500
		dims = [1, 2, 3]
		num_colours = [2, 3]
		for dim, s in product(dims, num_colours):
			points = rng.uniform(size=(dim, num_points))
			colours = rng.integers(0, s - 1, size=num_points).tolist()
			filtration, numerical_errors = ch.chromatic.delrips(points, colours, max_num_threads=0)
			assert not numerical_errors
			dgms_kchromatic = ch.sixpack.from_filtration(
				filtration,
				mapping_method=ch.sixpack.KChromaticInclusion(1),
			)
			dgms_kgluing = ch.sixpack.from_filtration(
				filtration,
				mapping_method=ch.sixpack.KChromaticGluingMap(1),
			)
			# For k=1, the k-chromatic inclusion and gluing map should yield the same result.
			dists = [
				bottleneck_distance(
					dgms_kchromatic.get_matrix(name, d),
					dgms_kgluing.get_matrix(name, d),
				)
				for name in get_args(ch.sixpack.DiagramEnsemble.DiagramName.__value__)
				for d in range(dim)
			]
			assert math.isclose(max(dists), 0, abs_tol=1e-6)
