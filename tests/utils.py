import math
from itertools import combinations

import chalc as ch


def assert_standard_simplex(
	n: int,
	filtration: ch.filtration.Filtration,
	*,
	is_filtered: bool,
	is_chromatic: bool,
) -> None:
	"""Assert that the filtration is a standard simplex of dimension n."""
	assert filtration.dimension == n
	assert len(filtration.simplices) == n + 1
	for i in range(n + 1):
		num_simplices = math.comb(n + 1, i + 1)
		assert list(filtration.simplices[i].keys()) == list(range(num_simplices)), (
			f"Dimension: {i}, Keys: {list(filtration.simplices[i].keys())}, "
			f"Wanted: {list(range(num_simplices))}"
		)
		gen_vertices = combinations(range(n + 1), i + 1)
		for k in range(num_simplices):
			assert filtration.simplices[i][k].dimension == i
			assert filtration.simplices[i][k].vertices == list(next(gen_vertices))
			assert filtration.simplices[i][k].label == k
			if not is_chromatic:
				assert filtration.simplices[i][k].colours == [0]
			if not is_filtered:
				assert filtration.simplices[i][k].filtration_value == 0
