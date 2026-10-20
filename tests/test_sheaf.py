"""Induced persistence maps use consistent interval bases and compose over F2."""

from collections.abc import Mapping
from itertools import combinations, product

import numpy as np
import pytest

from chalc.chromatic import delaunay_cech
from chalc.filtration import Filtration, standard_simplex
from chalc.morphisms import (
	AugmentedPersistenceModule,
	IndexArray,
	SixPack,
	SixPackMorphisms,
	SixPackPairings,
	SixPackRepresentatives,
	SubChromaticInclusion,
	sheaf,
)
from chalc.morphisms.sheaf import (
	ChromaticSheaf,
	PersistenceModuleMorphism,
	SubfiltrationInclusion,
	_induced_morphism,
	_solve,
)


def _rows(mapping: Mapping[int, list[int] | IndexArray]) -> dict[int, IndexArray]:
	return {key: np.asarray(value, dtype=np.int64) for key, value in mapping.items()}


def _comparison_filtrations(
	points: np.ndarray,
	colours: list[int],
) -> tuple[Filtration, Filtration]:
	return delaunay_cech(points, colours), delaunay_cech(points, [0] * points.shape[1])


def _point_inclusion(points: np.ndarray, colours: list[int]) -> PersistenceModuleMorphism:
	return _induced_morphism(*_comparison_filtrations(points, colours))


@pytest.mark.parametrize(
	("matrix", "rhs", "expected"),
	[
		({10: [2], 20: [4]}, [4, 2], [10, 20]),
		({10: [0, 2], 20: [1, 2]}, [1, 0], [10, 20]),
		({10: [0], 20: [2, 1], 30: [3, 0, 2]}, [1, 3], [10, 20, 30]),
		({10: [], 20: [1]}, [1], [20]),
		({10: [1], 20: [1], 30: []}, [1], [10]),
		({}, [], []),
		({10: [1]}, [], []),
	],
)
def test_solve(matrix: dict[int, list[int]], rhs: list[int], expected: list[int]) -> None:
	columns = _rows(matrix)
	vector = np.asarray(rhs, dtype=np.int64)
	original = {key: column.copy() for key, column in columns.items()}
	original_rhs = vector.copy()
	result = _solve(columns, vector)
	assert result.dtype == np.int64
	assert result.ndim == 1
	assert result.flags.c_contiguous
	np.testing.assert_equal(result, expected)
	np.testing.assert_equal(columns, original)
	np.testing.assert_equal(vector, original_rhs)


@pytest.mark.parametrize("matrix", [{}, {10: []}, {10: [0, 1], 20: [0, 1]}])
def test_solve_rejects_inconsistent_system(matrix: dict[int, list[int]]) -> None:
	with pytest.raises(ValueError, match="outside the column span"):
		_solve(_rows(matrix), np.asarray([0], dtype=np.int64))


def test_subfiltration_inclusion_allows_higher_dimensional_codomain() -> None:
	subfiltration = Filtration(3, 1)
	subfiltration.add_simplex([0, 1], 0)
	filtration = standard_simplex(2)
	inclusion = SubfiltrationInclusion(subfiltration, filtration)
	membership = [inclusion.simplex_in_domain(column) for column in filtration.boundary_matrix()]
	assert membership == [True, True, True, True, False, False, False]
	pack = inclusion.sixpack(get_representatives=True)
	np.testing.assert_array_equal(pack.domain.intervals(0), [[0, np.inf], [0, np.inf]])


@pytest.mark.parametrize(("source_vertices", "target_vertices"), [(3, 4), (4, 3)])
def test_subfiltration_inclusion_requires_same_vertex_universe(
	source_vertices: int,
	target_vertices: int,
) -> None:
	subfiltration = Filtration(source_vertices, 1)
	filtration = Filtration(target_vertices, 1)
	with pytest.raises(ValueError, match="same vertex universe"):
		SubfiltrationInclusion(subfiltration, filtration)


def _alive(module: AugmentedPersistenceModule, generator: int | np.integer, time: float) -> bool:
	birth, death = module.interval(generator)
	return birth <= time < death


def _square() -> np.ndarray:
	return np.array([[0.0, 1.0, 1.0, 0.0], [0.0, 0.0, 1.0, 1.0]])


@pytest.mark.parametrize("points", [_square(), np.array([[0.0, 1.0, 3.0]]), np.zeros((2, 1))])
def test_identity_inclusion(points: np.ndarray) -> None:
	bichromatic, monochromatic = _comparison_filtrations(points, [0] * points.shape[1])
	pack = SubChromaticInclusion(bichromatic, {0}).sixpack(get_representatives=True)
	assert isinstance(pack.representatives, SixPackRepresentatives)
	result = _induced_morphism(bichromatic, monochromatic)
	assert result.domain == result.codomain
	np.testing.assert_equal(
		result.map,
		{generator: [generator] for generator in pack.representatives.domain},
	)


def test_empty_source() -> None:
	result = ChromaticSheaf.from_chromatic_points(_square(), [1] * 4, 1).morphism(
		frozenset(),
		frozenset({1}),
	)
	assert not result.domain
	assert result.codomain
	np.testing.assert_equal(result.map, {})


def test_two_components_map_through_an_added_point() -> None:
	result = _point_inclusion(np.array([[0.0, 1.0, 3.0]]), [0, 1, 0])
	np.testing.assert_array_equal(result.domain.intervals(0), [[0, np.inf], [0, 1.5]])
	np.testing.assert_array_equal(result.codomain.intervals(0), [[0, np.inf], [0, 0.5], [0, 1.0]])
	np.testing.assert_equal(result.map, {0: [0], 1: [1, 2]})


def test_adding_centre_kills_the_square_loop() -> None:
	points = np.column_stack((_square(), [0.5, 0.5]))
	bichromatic, monochromatic = _comparison_filtrations(points, [0, 0, 0, 0, 1])
	pack = SubChromaticInclusion(bichromatic, {0}).sixpack(get_representatives=True)
	assert isinstance(pack.representatives, SixPackRepresentatives)
	result = _induced_morphism(bichromatic, monochromatic)
	loops = [
		generator
		for generator in pack.representatives.domain
		if result.domain.interval(generator)[0] == 0.5
	]
	assert len(loops) == 1
	assert result.domain.interval(loops[0]) == pytest.approx((0.5, np.sqrt(0.5)))
	assert loops[0] not in result.map
	assert result.codomain.intervals(1).shape == (0, 2)


def test_adding_distant_point_preserves_the_square_loop() -> None:
	points = np.column_stack((_square(), [3.0, 3.0]))
	bichromatic, monochromatic = _comparison_filtrations(points, [0, 0, 0, 0, 1])
	pack = SubChromaticInclusion(bichromatic, {0}).sixpack(get_representatives=True)
	assert isinstance(pack.representatives, SixPackRepresentatives)
	result = _induced_morphism(bichromatic, monochromatic)
	loops = [
		generator
		for generator in pack.representatives.domain
		if result.domain.interval(generator)[0] == 0.5
	]
	assert len(loops) == 1
	np.testing.assert_equal(len(result.map[loops[0]]), 1)
	target = result.map[loops[0]][0]
	assert result.codomain.interval(target) == result.domain.interval(loops[0])


def test_induced_morphism_uses_supplied_filtrations(monkeypatch: pytest.MonkeyPatch) -> None:
	points = np.random.default_rng(0).uniform(size=(2, 8))
	bichromatic, monochromatic = _comparison_filtrations(points, [0, 1] * 4)
	expected = _induced_morphism(bichromatic, monochromatic)

	def unexpected_construction(*_args: object, **_kwargs: object) -> Filtration:
		pytest.fail("The supplied filtrations should be reused.")

	monkeypatch.setattr(sheaf, "delaunay_cech", unexpected_construction)
	result = _induced_morphism(bichromatic, monochromatic)
	assert result == expected


def _rank_mod2(matrix: np.ndarray) -> int:
	matrix = matrix.copy()
	rank = 0
	for column in range(matrix.shape[1]):
		pivots = np.flatnonzero(matrix[rank:, column])
		if not len(pivots):
			continue
		pivot = rank + pivots[0]
		matrix[[rank, pivot]] = matrix[[pivot, rank]]
		for row in range(rank + 1, matrix.shape[0]):
			if matrix[row, column]:
				matrix[row] ^= matrix[rank]
		rank += 1
	return rank


def _assert_compatible_intervals(result: PersistenceModuleMorphism) -> None:
	for source, targets in result.map.items():
		birth, death = result.domain.interval(source)
		assert birth < death
		for target in targets:
			target_birth, target_death = result.codomain.interval(target)
			assert target_birth <= birth < target_death <= death


@pytest.mark.parametrize(("dimension", "seed"), [(1, 0), (2, 0), (2, 1), (2, 44), (3, 2)])
def test_induced_map_agrees_with_inclusion_of_cycles(dimension: int, seed: int) -> None:
	points = np.random.default_rng(seed).uniform(size=(dimension, 10))
	colours = [0, 1] * 5
	bichromatic, monochromatic = _comparison_filtrations(points, colours)
	result = _induced_morphism(bichromatic, monochromatic)
	pack_i = SubChromaticInclusion(bichromatic, {0}).sixpack(get_representatives=True)
	pack_j = SubfiltrationInclusion(monochromatic, bichromatic).sixpack(get_representatives=True)
	assert pack_i.codomain == pack_j.codomain
	assert isinstance(pack_i.representatives, SixPackRepresentatives)
	assert isinstance(pack_j.representatives, SixPackRepresentatives)
	columns = bichromatic.boundary_matrix()
	boundary = np.zeros((len(columns), len(columns)), dtype=np.uint8)
	for index, (facets, _, _, _) in enumerate(columns):
		boundary[facets, index] = 1
	_assert_compatible_intervals(result)

	# At every critical time, the image differs from the included source cycle
	# by a boundary, and the output cycles form a homology basis.
	times = sorted({0.0, *(column[2] for column in columns)})
	for time in times:
		available = [index for index, column in enumerate(columns) if column[2] <= time]
		boundaries = boundary[:, available]
		rank = _rank_mod2(boundaries)
		codomain_cycles = []
		for target, support in pack_j.representatives.domain.items():
			if _alive(result.codomain, target, time):
				cycle = np.zeros(len(columns), dtype=np.uint8)
				cycle[support] = 1
				codomain_cycles.append(cycle)
		augmented = np.column_stack((boundaries, *codomain_cycles))
		assert _rank_mod2(augmented) == rank + len(codomain_cycles)
		for source, support in pack_i.representatives.domain.items():
			if not _alive(result.domain, source, time):
				continue
			difference = np.zeros(len(columns), dtype=np.uint8)
			difference[support] = 1
			for target in result.map.get(source, []):
				if _alive(result.codomain, target, time):
					difference[pack_j.representatives.domain[int(target)]] ^= 1
			assert not np.any(boundary @ difference % 2)
			assert _rank_mod2(np.column_stack((boundaries, difference))) == rank


def _compose(
	first: PersistenceModuleMorphism,
	second: PersistenceModuleMorphism,
) -> dict[int, list[int]]:
	assert first.codomain == second.domain
	result = {}
	for source, intermediate in first.map.items():
		birth = first.domain.interval(source)[0]
		image = set()
		for middle in intermediate:
			for target in second.map.get(int(middle), ()):
				if _alive(second.codomain, target, birth):
					image.symmetric_difference_update({int(target)})
		if image:
			result[source] = sorted(image)
	return result


@pytest.mark.parametrize(("dimension", "seed"), product([1, 2, 3], [0, 44]))
def test_nested_inclusions_have_consistent_bases_and_compose(dimension: int, seed: int) -> None:
	points = np.random.default_rng(seed).uniform(size=(dimension, 14))
	# Keep the inherited order, including when a subset is not a prefix.
	middle = np.array([0, 2, 3, 5, 6, 8, 10, 11])
	inner = middle[::2]
	middle_colours = np.ones(14, dtype=int)
	middle_colours[middle] = 0
	inner_colours = np.ones(14, dtype=int)
	inner_colours[inner] = 0
	first = _point_inclusion(points[:, middle], [0, 1] * 4)
	second = _point_inclusion(points, middle_colours.tolist())
	direct = _point_inclusion(points, inner_colours.tolist())
	assert first.domain == direct.domain
	assert second.codomain == direct.codomain
	np.testing.assert_equal(_compose(first, second), direct.map)
	composed = second @ first
	assert composed == direct
	for time in _event_times(first.domain, first.codomain, second.codomain):
		np.testing.assert_array_equal(_matrix_at(composed, time), _matrix_at(direct, time))


def test_global_solve_discards_disjoint_interval_coefficients(
	monkeypatch: pytest.MonkeyPatch,
) -> None:
	# Nonzero intervals [0, 2), [1, 3), [2, 4) form the common barcode of j.
	pairings: dict[int, int | None] = {0: 2, 1: 4, 3: 5}
	times = [0.0, 1.0, 2.0, 2.0, 3.0, 4.0]

	def maps(mapping: dict[int, list[int]]) -> SixPackMorphisms:
		return SixPackMorphisms(
			kernel_to_domain=_rows({}),
			domain_to_image=_rows({}),
			image_to_codomain=_rows({}),
			domain_to_codomain=_rows(mapping),
			codomain_to_cokernel=_rows({}),
		)

	pack_i = SixPack(
		pairings=SixPackPairings(domain={3: 5}, codomain=pairings),
		entrance_times=times,
		degrees=[0] * 6,
		morphisms=maps({3: [3]}),
	)
	matrix_j = {0: [0], 1: [0, 1], 3: [1, 3]}
	pack_j = SixPack(
		pairings=SixPackPairings(domain=pairings, codomain=pairings),
		entrance_times=times,
		degrees=[0] * 6,
		morphisms=maps(matrix_j),
	)
	np.testing.assert_equal(_solve(_rows(matrix_j), np.asarray([3], dtype=np.int64)), [0, 1, 3])
	monkeypatch.setattr(SubChromaticInclusion, "sixpack", lambda _self, **_kwargs: pack_i)
	monkeypatch.setattr(SubfiltrationInclusion, "sixpack", lambda _self, **_kwargs: pack_j)
	filtration = delaunay_cech(np.zeros((1, 1)), [0])
	result = _induced_morphism(filtration, filtration)
	np.testing.assert_equal(result.map, {0: [1, 3]})


def _event_times(*modules: AugmentedPersistenceModule) -> list[float]:
	times = sorted(
		{0.0}
		| {
			float(time) for module in modules for time in module.entrance_times if np.isfinite(time)
		},
	)
	return [times[0] - 1, *times, times[-1] + 1]


def _matrix_at(morphism: PersistenceModuleMorphism, time: float) -> np.ndarray:
	sources = sorted(b for b in morphism.domain.pairings if _alive(morphism.domain, b, time))
	targets = sorted(b for b in morphism.codomain.pairings if _alive(morphism.codomain, b, time))
	rows = {birth: row for row, birth in enumerate(targets)}
	result = np.zeros((len(targets), len(sources)), dtype=np.uint8)
	for column, source in enumerate(sources):
		for target in morphism.map.get(source, []):
			if target in rows:
				result[rows[target], column] ^= 1
	return result


@pytest.fixture
def essential_module() -> AugmentedPersistenceModule:
	return AugmentedPersistenceModule({0: None, 1: None, 2: None}, [0.0] * 3, {0: 0, 1: 0, 2: 0})


@pytest.mark.parametrize("mapping", [{9: [0]}, {0: [9]}])
def test_morphism_rejects_invalid_generators(
	essential_module: AugmentedPersistenceModule,
	mapping: dict[int, list[int]],
) -> None:
	with pytest.raises(ValueError, match=r"refer to .* intervals"):
		PersistenceModuleMorphism(essential_module, essential_module, _rows(mapping))


def _interval_module(birth: float, death: float, degree: int = 0) -> AugmentedPersistenceModule:
	if np.isinf(death):
		return AugmentedPersistenceModule({0: None}, [birth], {0: degree})
	return AugmentedPersistenceModule({0: 1}, [birth, death], {0: degree})


@pytest.mark.parametrize(
	("source_interval", "target_interval"),
	[
		((1.0, 4.0), (0.0, 3.0)),
		((1.0, 4.0), (1.0, 4.0)),
		((1.0, np.inf), (0.0, np.inf)),
		((1.0, np.inf), (0.0, 2.0)),
		((1.0, 4.0), (0.0, 4.0)),
		((1.0, 4.0), (1.0, 2.0)),
	],
)
def test_morphism_accepts_compatible_interval_coefficients(
	source_interval: tuple[float, float],
	target_interval: tuple[float, float],
) -> None:
	source = _interval_module(*source_interval)
	target = _interval_module(*target_interval)
	result = PersistenceModuleMorphism(source, target, _rows({0: [0]}))
	np.testing.assert_array_equal(_matrix_at(result, source_interval[0]), [[1]])
	_assert_compatible_intervals(result)


@pytest.mark.parametrize(
	("source_interval", "target_interval"),
	[
		((0.0, 2.0), (1.0, 2.0)),  # Target is born too late.
		((1.0, 2.0), (0.0, 3.0)),  # Target dies too late.
		((2.0, 4.0), (0.0, 2.0)),  # Half-open intervals meet but do not overlap.
		((2.0, 4.0), (0.0, 1.0)),
		((0.0, 2.0), (2.0, 3.0)),
		((1.0, 2.0), (0.0, np.inf)),
	],
)
def test_morphism_rejects_incompatible_interval_coefficients(
	source_interval: tuple[float, float],
	target_interval: tuple[float, float],
) -> None:
	source = _interval_module(*source_interval)
	target = _interval_module(*target_interval)
	with pytest.raises(
		ValueError, match="target birth <= source birth < target death <= source death"
	):
		PersistenceModuleMorphism(source, target, _rows({0: [0]}))


@pytest.mark.parametrize(("source_degree", "target_degree"), [(0, 1), (1, 0)])
def test_morphism_requires_matching_homological_degrees(
	source_degree: int,
	target_degree: int,
) -> None:
	source = _interval_module(0.0, 1.0, source_degree)
	target = _interval_module(0.0, 1.0, target_degree)
	with pytest.raises(ValueError, match="preserve homological degree"):
		PersistenceModuleMorphism(source, target, _rows({0: [0]}))


def test_zero_morphism_allows_unrelated_intervals_and_degrees() -> None:
	source = _interval_module(0.0, 1.0, 0)
	target = _interval_module(2.0, 3.0, 1)
	result = PersistenceModuleMorphism.zero(source, target)
	assert result.domain == source
	assert result.codomain == target
	assert result.map == {}


def test_morphism_composition_over_f2(essential_module: AugmentedPersistenceModule) -> None:
	first = PersistenceModuleMorphism(
		essential_module,
		essential_module,
		_rows({0: [0, 1], 1: [1, 2], 2: [2]}),
	)
	second = PersistenceModuleMorphism(
		essential_module,
		essential_module,
		_rows({0: [0, 2], 1: [1, 2]}),
	)
	original_first = {key: value.copy() for key, value in first.map.items()}
	original_second = {key: value.copy() for key, value in second.map.items()}
	composed = second @ first
	assert composed.domain == first.domain
	assert composed.codomain == second.codomain
	np.testing.assert_array_equal(_matrix_at(composed, 0), [[1, 0, 0], [1, 1, 0], [0, 1, 0]])
	np.testing.assert_equal(first.map, original_first)
	np.testing.assert_equal(second.map, original_second)


def test_morphism_identity_and_zero(essential_module: AugmentedPersistenceModule) -> None:
	module = essential_module
	morphism = PersistenceModuleMorphism(module, module, _rows({0: [0, 2], 2: [1]}))
	identity = PersistenceModuleMorphism.identity(module)
	zero = PersistenceModuleMorphism.zero(module, module)
	np.testing.assert_equal(identity.map, {birth: [birth] for birth in module.pairings})
	assert zero.domain == module
	assert zero.codomain == module
	np.testing.assert_equal(zero.map, {})
	for composed in (identity @ morphism, morphism @ identity):
		np.testing.assert_array_equal(_matrix_at(composed, 0), _matrix_at(morphism, 0))
	# In particular, zero @ identity must treat omitted columns as zero.
	for composed in (zero @ identity, identity @ zero, zero @ morphism, morphism @ zero):
		assert not _matrix_at(composed, 0).any()
	empty = PersistenceModuleMorphism.identity(AugmentedPersistenceModule())
	np.testing.assert_equal((empty @ empty).map, {})


def test_morphism_rejects_incompatible_composition(
	essential_module: AugmentedPersistenceModule,
) -> None:
	other = AugmentedPersistenceModule({0: None}, [1.0], {0: 0})
	with pytest.raises(RuntimeError, match="Incompatible morphisms"):
		PersistenceModuleMorphism.identity(other) @ PersistenceModuleMorphism.identity(
			essential_module,
		)


@pytest.mark.parametrize("target_death", [1.5, 2.0])
def test_composition_discards_disjoint_interval_coefficients(target_death: float) -> None:
	# Both arrows are nonzero, but the source and final target do not overlap.
	source = AugmentedPersistenceModule({0: 1}, [2.0, 4.0], {0: 0})
	middle = AugmentedPersistenceModule({0: 1}, [1.0, 3.0], {0: 0})
	target = AugmentedPersistenceModule({0: 1}, [0.0, target_death], {0: 0})
	first = PersistenceModuleMorphism(source, middle, _rows({0: [0]}))
	second = PersistenceModuleMorphism(middle, target, _rows({0: [0]}))
	composed = second @ first
	assert composed.map == {}
	assert composed == PersistenceModuleMorphism.zero(source, target)
	for time in _event_times(source, middle, target):
		actual = _matrix_at(composed, time)
		np.testing.assert_array_equal(
			actual,
			_matrix_at(second, time) @ _matrix_at(first, time) % 2,
		)
		assert not actual.any()


def test_composition_filters_only_disjoint_interval_coefficients() -> None:
	source = _interval_module(2.0, 4.0)
	middle = _interval_module(1.0, 3.0)
	target = AugmentedPersistenceModule(
		{0: 4, 1: 5, 2: 6, 3: 7},
		[0.0, 0.0, 0.0, 0.0, 1.5, 2.0, 2.5, 3.0],
		dict.fromkeys(range(4), 0),
	)
	first = PersistenceModuleMorphism(source, middle, _rows({0: [0]}))
	second = PersistenceModuleMorphism(middle, target, _rows({0: [0, 1, 2, 3]}))
	composed = second @ first
	np.testing.assert_equal(composed.map, {0: [2, 3]})
	assert composed.map[0].dtype == np.int64
	_assert_compatible_intervals(composed)
	for time in _event_times(source, middle, target):
		np.testing.assert_array_equal(
			_matrix_at(composed, time),
			_matrix_at(second, time) @ _matrix_at(first, time) % 2,
		)


@pytest.mark.parametrize(
	("colours", "height"),
	[([], 0), ([7], 1), ([2, 7, 9], 0), ([2, 7, 9], 1), ([2, 7, 9], 2), ([2, 7, 9], 10)],
)
def test_sheaf_construction(colours: list[int], height: int) -> None:
	points = np.arange(len(colours), dtype=float).reshape((1, -1))
	result = ChromaticSheaf.from_chromatic_points(points, colours, height)
	all_colours = frozenset(colours)
	assert result.colours == all_colours
	assert result.max_num_colours == min(height, len(all_colours))
	assert not result.stalk(frozenset())
	for size in range(result.max_num_colours + 1):
		for subset in combinations(sorted(all_colours), size):
			colour_set = frozenset(subset)
			module = result.stalk(colour_set)
			identity = result.morphism(colour_set, colour_set)
			assert identity.domain == module
			assert identity.codomain == module
			np.testing.assert_equal(identity.map, {birth: [birth] for birth in module.pairings})
			assert result.morphism(colour_set, colour_set) is identity
			from_empty = result.morphism(frozenset(), colour_set)
			assert from_empty.domain == result.stalk(frozenset())
			assert from_empty.codomain == module
			np.testing.assert_equal(from_empty.map, {})
			if colour_set:
				in_set = [c in colour_set for c in colours]
				expected = _point_inclusion(points[:, in_set], [0] * sum(in_set)).domain
				assert module == expected


@pytest.mark.parametrize(
	("colours", "height", "message"),
	[([0], 1, "Length of `colours`"), ([0, 1], -1, "`max_num_colours`")],
)
def test_sheaf_rejects_invalid_inputs(colours: list[int], height: int, message: str) -> None:
	with pytest.raises(ValueError, match=message):
		ChromaticSheaf.from_chromatic_points(np.array([[0.0, 1.0]]), colours, height)


@pytest.mark.parametrize(
	("source", "target", "message"),
	[
		(frozenset({0}), frozenset({1}), "source <= target"),
		(frozenset(), frozenset({9}), "source <= target"),
		(frozenset(), frozenset({0, 1, 2}), "at most"),
	],
)
def test_sheaf_rejects_invalid_morphism_requests(
	source: frozenset[int],
	target: frozenset[int],
	message: str,
) -> None:
	result = ChromaticSheaf.from_chromatic_points(np.array([[0.0, 1.0, 3.0]]), [0, 1, 2], 2)
	with pytest.raises(ValueError, match=message):
		result.morphism(source, target)


@pytest.mark.parametrize("prefetch", [False, True])
def test_sheaf_composes_three_edges_and_caches(*, prefetch: bool) -> None:
	points = np.array([[0.0, 1.0, 3.0, 7.0]])
	result = ChromaticSheaf.from_chromatic_points(points, [0, 1, 2, 3], 4)
	source, target = frozenset({0}), frozenset({0, 1, 2, 3})
	if prefetch:
		result.morphism(source, frozenset({0, 1, 2}))
	actual = result.morphism(source, target)
	expected = _point_inclusion(points, [0, 1, 1, 1])
	assert actual.domain == result.stalk(source) == expected.domain
	assert actual.codomain == result.stalk(target) == expected.codomain
	for time in _event_times(actual.domain, actual.codomain):
		np.testing.assert_array_equal(_matrix_at(actual, time), _matrix_at(expected, time))
	assert result.morphism(source, target) is actual


@pytest.mark.parametrize("dimension", [1, 2])
def test_sheaf_maps_agree_with_direct_inclusions(dimension: int) -> None:
	points = np.random.default_rng(44).uniform(size=(dimension, 9))
	colours = [2, 7, 11] * 3
	result = ChromaticSheaf.from_chromatic_points(points, colours, 3)
	all_colours = sorted(result.colours)
	for size in range(1, 4):
		for subset in combinations(all_colours, size):
			target = frozenset(subset)
			in_target = [c in target for c in colours]
			for source_size in range(1, size + 1):
				for source_subset in combinations(subset, source_size):
					source = frozenset(source_subset)
					bridge_colours = [int(c not in source) for c in colours if c in target]
					expected = _point_inclusion(points[:, in_target], bridge_colours)
					actual = result.morphism(source, target)
					assert actual.domain == result.stalk(source) == expected.domain
					assert actual.codomain == result.stalk(target) == expected.codomain
					for time in _event_times(actual.domain, actual.codomain):
						np.testing.assert_array_equal(
							_matrix_at(actual, time),
							_matrix_at(expected, time),
						)
