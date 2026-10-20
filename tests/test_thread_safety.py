"""Concurrent uses supported by the thread-safety policy."""

from __future__ import annotations

import copy
import subprocess
import sys
import sysconfig
from concurrent.futures import ThreadPoolExecutor
from itertools import combinations
from threading import Barrier, Lock
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pytest

import chalc as ch

if TYPE_CHECKING:
	from collections.abc import Callable
	from pathlib import Path

_WORKERS = 4
_REPEATS = 4
_GEOMETRIES = ("delaunay", "delaunay_rips", "alpha", "delaunay_cech")


def _run_concurrently[T](operation: Callable[[int], T]) -> list[T]:
	barrier = Barrier(_WORKERS)

	def run(worker: int) -> T:
		barrier.wait(timeout=30)
		return operation(worker)

	with ThreadPoolExecutor(max_workers=_WORKERS) as pool:
		return list(pool.map(run, range(_WORKERS)))


@pytest.fixture(scope="module", autouse=True)
def _check_gil() -> None:
	if sysconfig.get_config_var("Py_GIL_DISABLED"):
		assert not sys._is_gil_enabled(), (  # noqa: SLF001
			"Thread-safety tests must run with the GIL disabled."
		)


def test_import_keeps_gil_disabled() -> None:
	if not sysconfig.get_config_var("Py_GIL_DISABLED"):
		pytest.skip("Requires a free-threaded Python build.")
	# Ignore PYTHON_GIL and PYTHONPATH so neither can conceal an incompatible import.
	result = subprocess.run(
		[
			sys.executable,
			"-E",
			"-c",
			"import chalc; import sys; assert not sys._is_gil_enabled()",
		],
		capture_output=True,
		text=True,
		check=False,
		timeout=60,
	)
	assert result.returncode == 0, result.stdout + result.stderr


def _points(seed: int) -> tuple[np.ndarray, np.ndarray]:
	rng = np.random.default_rng(seed)
	return rng.random((2, 24)), np.arange(24, dtype=np.uint16) % 2


def _geometry(
	method: str,
	points: np.ndarray,
	colours: np.ndarray | list[int],
	*,
	parallel: bool,
) -> ch.filtration.Filtration:
	if method == "delaunay":
		return ch.chromatic.delaunay(points, colours, parallel=parallel)
	return getattr(ch.chromatic, method)(points, colours, max_num_threads=2 if parallel else 1)


@pytest.mark.parametrize("method", _GEOMETRIES)
@pytest.mark.parametrize("parallel", [False, True])
@pytest.mark.parametrize("array_colours", [False, True])
def test_independent_geometry_matches_serial(
	method: str,
	*,
	parallel: bool,
	array_colours: bool,
) -> None:
	inputs = [_points(44 + worker) for worker in range(_WORKERS)]
	originals = [(points.copy(), colours.copy()) for points, colours in inputs]
	expected = [
		_geometry(method, points, colours, parallel=False).boundary_matrix()
		for points, colours in inputs
	]

	def compute(worker: int) -> None:
		points, colours = inputs[worker]
		colour_input = colours if array_colours else colours.tolist()
		for _ in range(_REPEATS):
			filtration = _geometry(method, points, colour_input, parallel=parallel)
			assert filtration.is_filtration()
			assert filtration.boundary_matrix() == expected[worker]

	_run_concurrently(compute)
	for (points, colours), (original_points, original_colours) in zip(
		inputs,
		originals,
		strict=True,
	):
		np.testing.assert_array_equal(points, original_points)
		np.testing.assert_array_equal(colours, original_colours)


@pytest.mark.parametrize("method", _GEOMETRIES)
def test_geometry_shares_unchanged_arrays(method: str) -> None:
	points, colours = _points(44)
	points.setflags(write=False)
	colours.setflags(write=False)
	expected = _geometry(method, points, colours, parallel=False).boundary_matrix()

	def compute(_worker: int) -> None:
		for _ in range(_REPEATS):
			assert _geometry(method, points, colours, parallel=True).boundary_matrix() == expected

	_run_concurrently(compute)


def _make_filtration(worker: int = 0) -> ch.filtration.Filtration:
	filtration = ch.filtration.complete_complex(7, 2)
	for label, vertex in filtration.simplices[0].items():
		vertex.filtration_value = float(label + worker)
		vertex.set_colour((label + worker) % 2)
	filtration.propagate_filt_values(0)
	filtration.propagate_colours()
	return filtration


def test_independent_filtration_updates() -> None:
	expected = [_make_filtration(worker).boundary_matrix() for worker in range(_WORKERS)]

	def compute(worker: int) -> None:
		for _ in range(_REPEATS):
			assert _make_filtration(worker).boundary_matrix() == expected[worker]

	_run_concurrently(compute)


def _filtration_snapshot(filtration: ch.filtration.Filtration) -> tuple[object, ...]:
	return (
		len(filtration),
		filtration.dimension,
		filtration.is_filtration(),
		[0, 1] in filtration,
		filtration.boundary_matrix(),
		filtration.skeleton(1).boundary_matrix(),
		tuple(tuple(sorted(simplices)) for simplices in filtration.simplices),
		sorted(
			(
				simplex.dimension,
				simplex.label,
				tuple(simplex.vertices),
				simplex.filtration_value,
				tuple(simplex.colours),
				tuple(facet.label for facet in simplex.facets),
			)
			for simplex in filtration
		),
	)


def test_shared_filtration_reads() -> None:
	filtration = _make_filtration()
	expected = _filtration_snapshot(filtration)

	def read(_worker: int) -> None:
		for _ in range(12):
			assert _filtration_snapshot(filtration) == expected

	_run_concurrently(read)
	assert _filtration_snapshot(filtration) == expected


def test_shared_filtration_updates_with_caller_lock() -> None:
	filtration = ch.filtration.Filtration(n=6, k=2)
	vertex = filtration.simplices[0][0]
	lock = Lock()
	triangles = list(combinations(range(6), 3))

	def use_filtration(worker: int) -> None:
		for iteration in range(20):
			with lock:
				if worker < 2:
					filtration.add_simplex(
						list(triangles[(worker + 2 * iteration) % len(triangles)]),
						vertex.filtration_value + 1,
					)
					vertex.filtration_value += 1
					vertex.set_colour(int(vertex.filtration_value) % 3)
					filtration.propagate_filt_values(0)
					filtration.propagate_colours()
				assert filtration.is_filtration()
				vertices = filtration.simplices[0]
				for simplex in filtration:
					assert simplex.filtration_value == max(
						vertices[label].filtration_value for label in simplex.vertices
					)
					assert set(simplex.colours) == {
						vertices[label].colours[0] for label in simplex.vertices
					}
				assert len(filtration.boundary_matrix()) == len(filtration)

	_run_concurrently(use_filtration)
	assert vertex.filtration_value == 40
	assert len(filtration) == 41


def _make_morphism(
	kind: str,
	filtration: ch.filtration.Filtration,
) -> ch.morphisms.FiltrationMorphism:
	if kind == "sub_inclusion":
		return ch.morphisms.SubChromaticInclusion(filtration, [{0}])
	if kind == "k_inclusion":
		return ch.morphisms.KChromaticInclusion(filtration, 1)
	if kind == "sub_gluing":
		return ch.morphisms.SubChromaticGluing(filtration, [[{0}], [{1}]])
	return ch.morphisms.KChromaticGluing(filtration, 1)


def _assert_same_pack(actual: ch.morphisms.SixPack, expected: ch.morphisms.SixPack) -> None:
	assert actual == expected
	assert actual.pairings == expected.pairings
	assert actual.representatives == expected.representatives
	assert actual.morphisms == expected.morphisms
	assert actual.max_nonempty_degree() == expected.max_nonempty_degree()
	for name, module in actual:
		reference = getattr(expected, name)
		assert hash(module) == hash(reference)
		np.testing.assert_array_equal(module.nonempty_degrees(), reference.nonempty_degrees())
		for degree, intervals in module.intervals().items():
			np.testing.assert_array_equal(intervals, reference.intervals(degree))


@pytest.mark.parametrize("kind", ["sub_inclusion", "k_inclusion", "sub_gluing", "k_gluing"])
def test_shared_morphism_sixpack_calls(kind: str) -> None:
	filtration = _make_filtration()
	before = _filtration_snapshot(filtration)
	morphism = _make_morphism(kind, filtration)
	expected = morphism.sixpack(get_representatives=True)

	def compute(_worker: int) -> None:
		for _ in range(_REPEATS):
			_assert_same_pack(morphism.sixpack(get_representatives=True), expected)
			other_morphism = _make_morphism(kind, filtration)
			_assert_same_pack(other_morphism.sixpack(get_representatives=True), expected)

	_run_concurrently(compute)
	assert _filtration_snapshot(filtration) == before


def test_shared_sixpack_reads() -> None:
	pack = _make_morphism("sub_inclusion", _make_filtration()).sixpack(get_representatives=True)
	expected = copy.deepcopy(pack)

	def read(_worker: int) -> None:
		for _ in range(12):
			_assert_same_pack(pack, expected)

	_run_concurrently(read)
	_assert_same_pack(pack, expected)


def test_sixpack_hdf5_separate_files(tmp_path: Path) -> None:
	pack = _make_morphism("sub_inclusion", _make_filtration()).sixpack(get_representatives=True)

	def save_and_load(worker: int) -> None:
		path = tmp_path / f"sixpack-{worker}.h5"
		for _ in range(_REPEATS):
			with h5py.File(path, "w") as file:
				pack.to_hdf5(file)
			with h5py.File(path, "r") as file:
				loaded = ch.morphisms.SixPack.from_hdf5(file)
			_assert_same_pack(loaded, pack)

	_run_concurrently(save_and_load)
