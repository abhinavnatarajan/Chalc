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
