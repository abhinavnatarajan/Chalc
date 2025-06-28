import math
from itertools import combinations

import chalc as ch


def assert_clique_complex(
	filtration: ch.filtration.Filtration,
	n: int,
	k: int,
	*,
	is_filtered: bool,
	is_chromatic: bool,
) -> None:
	"""Assert that the filtration is the k-skeleton of the standard simplex of dimension n."""
	assert filtration.dimension == k
	assert len(filtration.simplices) == k + 1
	for i in range(k + 1):
		num_simplices = math.comb(n + 1, i + 1)
		assert sorted(filtration.simplices[i].keys()) == list(range(num_simplices)), (
			f"Dimension: {i}, Keys: {list(filtration.simplices[i].keys())}, "
			f"Wanted: {list(range(num_simplices))}"
		)
		gen_vertices = combinations(range(n + 1), i + 1)
		for j in range(num_simplices):
			assert filtration.simplices[i][j].dimension == i
			assert sorted(filtration.simplices[i][j].vertices) == list(next(gen_vertices))
			assert filtration.simplices[i][j].label == j
			if not is_chromatic:
				assert filtration.simplices[i][j].colours == [0]
			if not is_filtered:
				assert filtration.simplices[i][j].filtration_value == 0.0
