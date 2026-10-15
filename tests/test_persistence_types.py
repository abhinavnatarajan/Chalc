"""Regression tests for augmented persistence modules and 6-packs."""

import gc
import weakref
from collections.abc import Mapping
from dataclasses import fields
from pathlib import Path

import h5py
import numpy as np
import pytest

from chalc.morphisms import (
	AugmentedPersistenceModule,
	IndexArray,
	SixPack,
	SixPackMorphisms,
	SixPackPairings,
	SixPackRepresentatives,
)
from chalc.morphisms.types import _endpoint_indices


def _rows(mapping: Mapping[int, list[int] | IndexArray]) -> dict[int, IndexArray]:
	return {key: np.asarray(value, dtype=np.int64) for key, value in mapping.items()}


@pytest.fixture
def module() -> AugmentedPersistenceModule:
	return AugmentedPersistenceModule(
		{0: 2, 1: None, 3: 4},
		[0.0, 0.0, 1.0, 2.0, 3.0],
		{0: 0, 1: 2, 3: 2},
	)


def test_intervals_by_degree(module: AugmentedPersistenceModule) -> None:
	np.testing.assert_equal(module.nonempty_degrees(), [0, 2])
	assert module.nonempty_degrees().dtype == np.int64
	assert module.nonempty_degrees().flags.c_contiguous
	assert list(module.intervals(np.array([2, 0], dtype=np.int64))) == [2, 0]
	np.testing.assert_array_equal(module.entrance_times, [0.0, 0.0, 1.0, 2.0, 3.0])
	assert module.entrance_times.dtype == np.float64
	np.testing.assert_array_equal(module.intervals(0), [[0.0, 1.0]])
	np.testing.assert_array_equal(module.intervals(2), [[0.0, np.inf], [2.0, 3.0]])
	assert module.intervals(1).shape == (0, 2)
	assert module.intervals(1).dtype == np.float64
	assert list(module.intervals()) == [0, 2]
	assert list(module.intervals([2, 1, 0])) == [2, 1, 0]
	np.testing.assert_array_equal(module.intervals(np.int64(2)), module.intervals(2))


def test_empty_module() -> None:
	module = AugmentedPersistenceModule()
	assert not module
	np.testing.assert_equal(module.nonempty_degrees(), [])
	assert module.nonempty_degrees().dtype == np.int64
	assert module.nonempty_degrees().shape == (0,)
	assert module.entrance_times.shape == (0,)
	assert module.entrance_times.dtype == np.float64
	assert module.intervals() == {}
	assert module.intervals(0).shape == (0, 2)


@pytest.mark.parametrize(
	("birth_index", "expected"),
	[(0, (0.0, 1.0)), (1, (0.0, np.inf)), (3, (2.0, 3.0))],
)
def test_single_interval(
	module: AugmentedPersistenceModule,
	birth_index: int,
	expected: tuple[float, float],
) -> None:
	assert module.interval(birth_index) == expected


@pytest.mark.parametrize("birth_index", [-1, 2, 4, 5])
def test_single_interval_requires_generator(
	module: AugmentedPersistenceModule,
	birth_index: int,
) -> None:
	with pytest.raises(KeyError):
		module.interval(birth_index)


@pytest.mark.parametrize("degrees", ["0", [0, "1"], 0.5, {0}])
def test_invalid_requested_degrees(module: AugmentedPersistenceModule, degrees: object) -> None:
	with pytest.raises(TypeError, match="degrees"):
		# pyrefly: ignore [no-matching-overload]
		module.intervals(degrees)


@pytest.mark.parametrize(
	("intervals", "times", "degrees"),
	[
		({0: 1}, [1.0, 0.0], {0: 0}),
		({0: None}, [np.nan], {0: 0}),
		({0: 2}, [0.0, 1.0], {0: 0}),
		({-1: None}, [0.0], {-1: 0}),
		({1: 0}, [0.0, 1.0], {1: 0}),
		({0: None}, [0.0], {}),
		({0: None}, [0.0], {0: -1}),
	],
)
def test_invalid_module(
	intervals: dict[int, int | None],
	times: list[float],
	degrees: dict[int, int],
) -> None:
	with pytest.raises(ValueError, match=r"Endpoint|Interval|Degrees"):
		AugmentedPersistenceModule(intervals, times, degrees)


def test_module_equality_and_input_copies(module: AugmentedPersistenceModule) -> None:
	times = np.array([0.0, 1.0])
	intervals = {0: 1}
	degrees = {0: 0}
	original = AugmentedPersistenceModule(intervals, times, degrees)
	times[:] = 7.0
	intervals.clear()
	degrees.clear()
	equal = AugmentedPersistenceModule({0: 1}, [0.0, 1.0], {0: 0})
	assert original == equal
	assert hash(original) == hash(equal)
	assert original != module
	assert original != object()
	assert original != AugmentedPersistenceModule({0: None}, [0.0], {0: 0})


@pytest.mark.parametrize(
	"original",
	[
		AugmentedPersistenceModule(),
		AugmentedPersistenceModule(
			{3: 4, 0: 2, 1: None},
			[0.0, 0.0, 1.0, 2.0, 3.0],
			{1: 2, 3: 2, 0: 0},
		),
		AugmentedPersistenceModule({2: 2}, [0.0, 0.0, 1.0], {2: 1}),
	],
)
def test_module_hdf5_roundtrip(tmp_path: Path, original: AugmentedPersistenceModule) -> None:
	with h5py.File(tmp_path / "module.h5", "w") as file:
		original.to_hdf5(file.create_group("module"))
	with h5py.File(tmp_path / "module.h5", "r") as file:
		group = file["module"]
		assert group["pairings"].ndim == 2
		assert group["pairings"].shape[1] == 2
		assert group["pairings"].dtype == np.int64
		assert group["entrance_times"].dtype == np.float64
		assert group["degrees"].dtype == np.int64
		rows = group["pairings"][:]
		degrees = group["degrees"][:]
		loaded = AugmentedPersistenceModule(
			{int(b): int(d) if d != -1 else None for b, d in rows},
			group["entrance_times"][:],
			{int(b): int(degree) for (b, _), degree in zip(rows, degrees, strict=True)},
		)
	assert loaded == original


@pytest.fixture
def sixpack() -> SixPack:
	return SixPack(
		pairings=SixPackPairings(domain={1: 4, 2: None}, kernel={3: 5}, relative={3: None}),
		entrance_times=[0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0],
		degrees=[0, 0, 0, 1, 1, 2, 3],
	)


def test_sixpack_modules(sixpack: SixPack) -> None:
	assert set(SixPack.module_names()) == {
		"kernel",
		"cokernel",
		"domain",
		"codomain",
		"image",
		"relative",
	}
	assert [name for name, _ in sixpack] == list(SixPack.module_names())
	assert len(list(sixpack)) == 6
	assert all(module is getattr(sixpack, name) for name, module in sixpack)
	assert sixpack.domain is sixpack.domain
	assert sixpack
	assert not SixPack()
	assert SixPack().max_nonempty_degree() == -1
	assert sixpack.max_nonempty_degree() == 1
	assert sixpack.pairings.domain == {1: 4, 2: None}
	np.testing.assert_equal(sixpack.domain.nonempty_degrees(), [0])
	np.testing.assert_array_equal(sixpack.domain.intervals(0), [[0.0, 2.0], [0.0, np.inf]])
	np.testing.assert_equal(sixpack.kernel.nonempty_degrees(), [0])
	np.testing.assert_array_equal(sixpack.kernel.intervals(0), [[1.0, 3.0]])
	assert not sixpack.image
	assert sixpack.representatives is None
	assert sixpack.morphisms is None
	assert sixpack.metadata is None
	with pytest.raises(TypeError):
		hash(sixpack)


def test_sixpack_equality(sixpack: SixPack) -> None:
	assert sixpack != object()
	assert SixPack(entrance_times=[0.0], degrees=[0]) == SixPack()
	first = SixPack(pairings=SixPackPairings(domain={0: None}), entrance_times=[0.0], degrees=[0])
	assert first != SixPack()
	assert first != SixPack(
		pairings=SixPackPairings(domain={0: None}),
		entrance_times=[1.0],
		degrees=[0],
	)
	assert first != SixPack(
		pairings=SixPackPairings(domain={0: None}),
		entrance_times=[0.0],
		degrees=[1],
	)


def test_sixpack_reindexes_shared_endpoint_once() -> None:
	pack = SixPack(
		pairings=SixPackPairings(domain={1: 2, 2: 3}),
		entrance_times=[0.0, 0.0, 1.0, 2.0],
		degrees=[0, 0, 0, 0],
	)
	np.testing.assert_array_equal(pack.domain.entrance_times, [0.0, 1.0, 2.0])
	assert dict(pack.domain.pairings) == {0: 1, 1: 2}
	assert pack.domain.interval(0) == (0.0, 1.0)
	assert pack.domain.interval(1) == (1.0, 2.0)


def test_sixpack_releases_shared_endpoint_arrays() -> None:
	times = np.arange(100, dtype=np.float64)
	degrees = np.zeros(100, dtype=np.int64)
	times_ref = weakref.ref(times)
	degrees_ref = weakref.ref(degrees)
	pack = SixPack(
		pairings=SixPackPairings(domain={3: 7}),
		entrance_times=times,
		degrees=degrees,
	)
	del times, degrees
	gc.collect()
	assert times_ref() is None
	assert degrees_ref() is None
	assert pack.domain.interval(0) == (3.0, 7.0)


@pytest.mark.parametrize("pairings", [{-1: None}, {2: None}, {1: 0}, {0: 2}, {0: -1}])
def test_sixpack_rejects_invalid_raw_indices(pairings: dict[int, int | None]) -> None:
	with pytest.raises(ValueError, match="Invalid pairings"):
		SixPack(
			pairings=SixPackPairings(domain=pairings),
			entrance_times=[0.0, 1.0],
			degrees=[0, 0],
		)


def test_sixpack_takes_ownership_and_reindexes_representatives_and_maps() -> None:
	pairings = SixPackPairings(
		kernel={5: None},
		domain={2: None, 4: None},
		image={4: None},
		codomain={3: None},
	)
	representatives = SixPackRepresentatives(
		kernel=_rows({5: [10, 11]}),
		domain=_rows({2: [10], 4: [11]}),
		image=_rows({4: [11]}),
		codomain=_rows({3: [12]}),
		cokernel=_rows({}),
		relative=_rows({}),
	)
	morphisms = SixPackMorphisms(
		kernel_to_domain=_rows({5: [2, 4]}),
		domain_to_image=_rows({2: [], 4: [4]}),
		image_to_codomain=_rows({4: [3]}),
		domain_to_codomain=_rows({2: [], 4: [3]}),
		codomain_to_cokernel=_rows({3: []}),
	)

	kernel_support = representatives.kernel[5]
	domain_support = representatives.domain[2]
	pack = SixPack(
		pairings=pairings,
		representatives=representatives,
		morphisms=morphisms,
		entrance_times=[0.0, 0.0, 1.0, 1.0, 2.0, 3.0],
		degrees=[0, 0, 0, 0, 0, 1],
	)
	assert pack.pairings is pairings
	assert pack.representatives is representatives
	assert pack.morphisms is morphisms
	np.testing.assert_equal(representatives.kernel, {0: [10, 11]})
	np.testing.assert_equal(representatives.domain, {0: [10], 1: [11]})
	assert representatives.kernel[0] is kernel_support
	assert representatives.domain[0] is domain_support
	np.testing.assert_equal(morphisms.kernel_to_domain, {0: [0, 1]})
	np.testing.assert_equal(morphisms.domain_to_image, {0: [], 1: [0]})
	np.testing.assert_equal(morphisms.domain_to_codomain, {0: [], 1: [0]})


@pytest.mark.parametrize("empty", [False, True])
def test_hdf5_roundtrip(tmp_path: Path, sixpack: SixPack, *, empty: bool) -> None:
	original = SixPack() if empty else sixpack
	with h5py.File(tmp_path / "sixpack.h5", "w") as file:
		group = file.create_group("sixpack")
		original.to_hdf5(group)
		assert group["image/pairings"].shape == (0, 2)
		assert group["image/pairings"].dtype == np.int64
	with h5py.File(tmp_path / "sixpack.h5", "r") as file:
		loaded = SixPack.from_hdf5(file["sixpack"])
	assert loaded == original
	assert loaded.pairings == original.pairings
	assert loaded._endpoint_relabelling == original._endpoint_relabelling  # noqa: SLF001
	assert dict(loaded) == dict(original)
	assert loaded.representatives is None
	assert loaded.morphisms is None
	assert loaded.metadata is None


@pytest.mark.parametrize("damage", ["missing", "group", "dtype", "shape"])
def test_hdf5_invalid_dataset(tmp_path: Path, damage: str) -> None:
	with h5py.File(tmp_path / "invalid.h5", "w") as file:
		SixPack().to_hdf5(file)
		if damage == "missing":
			del file["domain/degrees"]
		elif damage == "group":
			del file["domain/pairings"]
			file.create_group("domain/pairings")
		elif damage == "dtype":
			del file["domain/pairings"]
			file["domain/pairings"] = np.empty((0, 2), dtype=np.float64)
		else:
			del file["domain/pairings"]
			file["domain/pairings"] = np.empty((0,), dtype=np.int64)
		with pytest.raises(RuntimeError, match="Invalid"):
			SixPack.from_hdf5(file)


@pytest.fixture
def sixpack_with_representatives() -> SixPack:
	return SixPack(
		pairings=SixPackPairings(
			kernel={5: None},
			domain={2: None, 4: None},
			image={4: None},
			codomain={3: None},
		),
		entrance_times=[0.0, 0.0, 1.0, 1.0, 2.0, 3.0],
		degrees=[0, 0, 0, 0, 0, 1],
		representatives=SixPackRepresentatives(
			kernel=_rows({5: [10, 11]}),
			domain=_rows({2: [10], 4: [11]}),
			image=_rows({4: [11]}),
			codomain=_rows({3: [12]}),
			cokernel=_rows({}),
			relative=_rows({}),
		),
		morphisms=SixPackMorphisms(
			kernel_to_domain=_rows({5: [2, 4]}),
			domain_to_image=_rows({2: [], 4: [4]}),
			image_to_codomain=_rows({4: [3]}),
			domain_to_codomain=_rows({2: [], 4: [3]}),
			codomain_to_cokernel=_rows({3: []}),
		),
	)


def test_hdf5_representatives_and_maps_roundtrip(
	tmp_path: Path,
	sixpack_with_representatives: SixPack,
) -> None:
	original = sixpack_with_representatives
	with h5py.File(tmp_path / "with_representatives.h5", "w") as file:
		original.to_hdf5(file)
	with h5py.File(tmp_path / "with_representatives.h5", "r") as file:
		loaded = SixPack.from_hdf5(file)
	assert loaded == original
	assert loaded.pairings == original.pairings
	assert loaded._endpoint_relabelling == original._endpoint_relabelling  # noqa: SLF001
	assert loaded.representatives == original.representatives
	assert loaded.morphisms == original.morphisms
	assert loaded.representatives is not None
	assert loaded.morphisms is not None
	np.testing.assert_equal(loaded.representatives.kernel, {0: [10, 11]})
	np.testing.assert_equal(loaded.morphisms.kernel_to_domain, {0: [0, 1]})
	np.testing.assert_equal(loaded.morphisms.domain_to_codomain, {0: [], 1: [0]})


def test_metadata_is_a_direct_consumer_attribute() -> None:
	pack = SixPack()
	metadata = {"value": object()}
	pack.metadata = metadata
	assert pack.metadata is metadata
	pack.metadata["extra"] = 7
	assert metadata["extra"] == 7


def test_hdf5_leaves_metadata_to_consumers(tmp_path: Path) -> None:
	pack = SixPack()
	pack.metadata = {"value": object()}
	with h5py.File(tmp_path / "metadata.h5", "w") as file:
		file.create_group("metadata").attrs["owner"] = "consumer"
		pack.to_hdf5(file)
		assert file["metadata"].attrs["owner"] == "consumer"
	with h5py.File(tmp_path / "metadata.h5", "r") as file:
		loaded = SixPack.from_hdf5(file)
	assert loaded.metadata is None


@pytest.mark.parametrize("with_representatives", [False, True])
@pytest.mark.parametrize("with_morphisms", [False, True])
def test_hdf5_optional_empty_representatives_and_maps(
	tmp_path: Path,
	*,
	with_representatives: bool,
	with_morphisms: bool,
) -> None:
	representatives = (
		SixPackRepresentatives(
			**{
				name: _rows(rows)
				for name, rows in ({name: {} for name in SixPack.module_names()}).items()
			},
		)
		if with_representatives
		else None
	)
	morphisms = (
		SixPackMorphisms(
			kernel_to_domain=_rows({}),
			domain_to_image=_rows({}),
			image_to_codomain=_rows({}),
			domain_to_codomain=_rows({}),
			codomain_to_cokernel=_rows({}),
		)
		if with_morphisms
		else None
	)
	with h5py.File(tmp_path / "optional_data.h5", "w") as file:
		SixPack(representatives=representatives, morphisms=morphisms).to_hdf5(file)
	with h5py.File(tmp_path / "optional_data.h5", "r") as file:
		loaded = SixPack.from_hdf5(file)
	assert loaded.representatives == representatives
	assert loaded.morphisms == morphisms


@pytest.mark.parametrize("damage", ["offsets", "keys", "indices", "group", "missing"])
def test_hdf5_rejects_invalid_sparse_data(
	tmp_path: Path,
	sixpack_with_representatives: SixPack,
	damage: str,
) -> None:
	with h5py.File(tmp_path / "invalid_sparse.h5", "w") as file:
		sixpack_with_representatives.to_hdf5(file)
		group = file["representatives/kernel"]
		if damage == "offsets":
			group["offsets"][-1] = 100
		elif damage == "keys":
			group["keys"][0] = 100
		elif damage == "indices":
			group["indices"][0] = -1
		elif damage == "group":
			del file["representatives"]
			file["representatives"] = np.array([], dtype=np.int64)
		else:
			del group["indices"]
		with pytest.raises(RuntimeError, match="Invalid"):
			SixPack.from_hdf5(file)


@pytest.mark.parametrize(
	"damage",
	["negative_birth", "death", "duplicate_birth", "times_length", "degrees_length"],
)
def test_hdf5_rejects_invalid_contents(tmp_path: Path, sixpack: SixPack, damage: str) -> None:
	with h5py.File(tmp_path / "invalid_contents.h5", "w") as file:
		sixpack.to_hdf5(file)
		if damage == "negative_birth":
			file["domain/pairings"][0, 0] = -1
		elif damage == "death":
			file["domain/pairings"][0, 1] = -2
		elif damage == "duplicate_birth":
			file["domain/pairings"][0, 0] = file["domain/pairings"][1, 0]
		elif damage == "times_length":
			del file["domain/entrance_times"]
			file["domain/entrance_times"] = np.array([], dtype=np.float64)
		else:
			del file["domain/degrees"]
			file["domain/degrees"] = np.array([], dtype=np.int64)
		with pytest.raises(RuntimeError, match="Invalid"):
			SixPack.from_hdf5(file)


@pytest.mark.parametrize("name", SixPack.module_names())
def test_sixpack_discards_degenerate_intervals_and_representatives(name: str) -> None:
	pairings = SixPackPairings(**{name: {0: 1, 1: 2, 3: None, 4: 4}})
	support = np.array([10, 11], dtype=np.int64)
	representatives = SixPackRepresentatives(
		**{
			name: _rows(rows)
			for name, rows in (
				{
					module: ({0: support, 1: [12], 3: [13], 4: [14]} if module == name else {})
					for module in SixPack.module_names()
				}
			).items()
		},
	)
	pack = SixPack(
		pairings=pairings,
		entrance_times=[0.0, 1.0, 1.0, 2.0, 3.0],
		degrees=[int(name == "kernel")] * 5,
		representatives=representatives,
	)
	module = getattr(pack, name)
	assert pack.pairings is pairings
	assert getattr(pairings, name) == {0: 1, 3: None}
	assert dict(module.pairings) == {0: 1, 2: None}
	np.testing.assert_array_equal(module.entrance_times, [0.0, 1.0, 2.0])
	np.testing.assert_equal(getattr(representatives, name), {0: [10, 11], 2: [13]})
	assert getattr(representatives, name)[0] is support
	np.testing.assert_equal(module.nonempty_degrees(), [0])
	assert pack._endpoint_relabelling == {  # noqa: SLF001
		module_name: ({0: 0, 3: 2} if module_name == name else {})
		for module_name in SixPack.module_names()
	}


@pytest.mark.parametrize("with_representatives", [False, True])
def test_sixpack_filters_morphism_sources_and_targets(
	tmp_path: Path,
	*,
	with_representatives: bool,
) -> None:
	pairings = {4: 5, 3: None, 1: 2, 0: 1}
	representatives = (
		SixPackRepresentatives(
			**{
				name: _rows(rows)
				for name, rows in (
					{
						name: (
							{birth: [10 + birth] for birth in pairings}
							if name in {"domain", "image", "codomain"}
							else {}
						)
						for name in SixPack.module_names()
					}
				).items()
			},
		)
		if with_representatives
		else None
	)
	morphisms = SixPackMorphisms(
		kernel_to_domain=_rows({}),
		domain_to_image=_rows({0: [1, 4], 1: [0], 3: [3], 4: [0]}),
		image_to_codomain=_rows({0: [0], 1: [1], 3: [3], 4: [4]}),
		domain_to_codomain=_rows({0: [0, 1, 4], 1: [3], 3: [1, 3], 4: [0]}),
		codomain_to_cokernel=_rows({0: [], 1: [], 3: [], 4: []}),
	)
	pack = SixPack(
		pairings=SixPackPairings(domain=pairings, image=pairings, codomain=pairings),
		entrance_times=[0.0, 1.0, 1.0, 2.0, 3.0, 3.0],
		degrees=[0] * 6,
		representatives=representatives,
		morphisms=morphisms,
	)
	assert pack.morphisms is morphisms
	np.testing.assert_equal(morphisms.domain_to_codomain, {0: [0], 2: [2]})
	np.testing.assert_equal(morphisms.domain_to_image, {0: [], 2: [2]})
	np.testing.assert_equal(morphisms.image_to_codomain, {0: [0], 2: [2]})
	np.testing.assert_equal(morphisms.codomain_to_cokernel, {0: [], 2: []})
	with h5py.File(tmp_path / "filtered_maps.h5", "w") as file:
		pack.to_hdf5(file)
	with h5py.File(tmp_path / "filtered_maps.h5", "r") as file:
		loaded = SixPack.from_hdf5(file)
	assert loaded == pack
	assert loaded.pairings == pack.pairings
	assert loaded._endpoint_relabelling == pack._endpoint_relabelling  # noqa: SLF001
	assert loaded.representatives == pack.representatives
	assert loaded.morphisms == pack.morphisms


@pytest.mark.parametrize("name", SixPack.module_names())
def test_sixpack_with_only_degenerate_intervals_is_empty(name: str) -> None:
	pairings: dict[int, int | None] = {0: 1, 2: 2}
	pack = SixPack(
		pairings=SixPackPairings(**{name: pairings}),
		entrance_times=[0.0, 0.0, 1.0],
		degrees=[int(name == "kernel")] * 3,
	)
	assert not pack
	assert pack.max_nonempty_degree() == -1
	assert getattr(pack.pairings, name) == {}
	assert getattr(pack, name).entrance_times.size == 0


@pytest.mark.parametrize("times", [[1.0, 0.0], [0.0, np.nan], [np.nan, 0.0]])
def test_sixpack_validates_times_before_filtering(times: list[float]) -> None:
	with pytest.raises(ValueError, match=r"nondecreasing.*NaN"):
		SixPack(pairings=SixPackPairings(domain={0: 1}), entrance_times=times, degrees=[0, 0])


@pytest.mark.parametrize(("name", "degree"), [("domain", -1), ("kernel", 0)])
def test_sixpack_validates_degrees_before_filtering(name: str, degree: int) -> None:
	pairings: dict[int, int | None] = {0: 1}
	with pytest.raises(ValueError, match="nonnegative"):
		SixPack(
			pairings=SixPackPairings(**{name: pairings}),
			entrance_times=[0.0, 0.0],
			degrees=[degree, degree],
		)


def test_hdf5_roundtrip_after_discarding_degenerate_intervals(tmp_path: Path) -> None:
	original = SixPack(
		pairings=SixPackPairings(domain={4: 5, 3: None, 1: 2, 0: 1}),
		entrance_times=[0.0, 1.0, 1.0, 2.0, 3.0, 3.0],
		degrees=[2, 0, 0, 1, 0, 0],
	)
	with h5py.File(tmp_path / "filtered.h5", "w") as file:
		original.to_hdf5(file)
		np.testing.assert_array_equal(file["domain/pairings"][:], [[3, -1], [0, 1]])
		np.testing.assert_array_equal(file["domain/entrance_times"][:], [0.0, 1.0, 2.0])
		np.testing.assert_array_equal(file["domain/degrees"][:], [1, 2])
	with h5py.File(tmp_path / "filtered.h5", "r") as file:
		loaded = SixPack.from_hdf5(file)
	assert loaded == original
	assert loaded.pairings == original.pairings
	assert loaded._endpoint_relabelling == original._endpoint_relabelling  # noqa: SLF001


@pytest.mark.parametrize("death", [1, None])
def test_hdf5_rejects_degenerate_intervals(tmp_path: Path, death: int | None) -> None:
	with h5py.File(tmp_path / "degenerate.h5", "w") as file:
		SixPack().to_hdf5(file)
		for name in ("pairings", "entrance_times", "degrees"):
			del file["domain"][name]
		times = [np.inf] if death is None else [1.0, 1.0]
		AugmentedPersistenceModule({0: death}, times, {0: 0}).to_hdf5(file["domain"])
		with pytest.raises(RuntimeError, match="Degenerate intervals"):
			SixPack.from_hdf5(file)


def test_sixpack_discards_interval_born_at_infinity() -> None:
	pack = SixPack(
		pairings=SixPackPairings(domain={0: None}),
		entrance_times=[np.inf],
		degrees=[0],
	)
	assert not pack
	assert pack.pairings.domain == {}


def test_hdf5_rebuilds_generator_relabelling_from_pairings(tmp_path: Path) -> None:
	with h5py.File(tmp_path / "relabel.h5", "w") as file:
		SixPack().to_hdf5(file)
		del file["domain"]
		# Pairing row order differs from endpoint order, and distinct endpoints
		# have equal times. Reconstruct by indices rather than by times.
		AugmentedPersistenceModule(
			{1: 3, 0: 2},
			[0.0, 0.0, 1.0, 1.0],
			{1: 1, 0: 0},
		).to_hdf5(file.create_group("domain"), pairings={11: 17, 5: 13})
	with h5py.File(tmp_path / "relabel.h5", "r") as file:
		loaded = SixPack.from_hdf5(file)
	assert loaded._endpoint_relabelling == {  # noqa: SLF001
		name: ({11: 1, 5: 0} if name == "domain" else {}) for name in SixPack.module_names()
	}
	assert dict(loaded.domain.pairings) == {1: 3, 0: 2}
	with h5py.File(tmp_path / "resaved.h5", "w") as file:
		loaded.to_hdf5(file)
	with h5py.File(tmp_path / "resaved.h5", "r") as file:
		reloaded = SixPack.from_hdf5(file)
	assert reloaded == loaded
	assert reloaded._endpoint_relabelling == loaded._endpoint_relabelling  # noqa: SLF001


@pytest.mark.parametrize(
	("pairings", "expected"),
	[({}, []), ({4: None, 1: 4, 2: 2}, [1, 2, 4]), ({9: 12, 3: 9}, [3, 9, 12])],
)
def test_endpoint_indices_are_sorted_unique_arrays(
	pairings: dict[int, int | None],
	expected: list[int],
) -> None:
	actual = _endpoint_indices(pairings)
	np.testing.assert_array_equal(actual, expected)
	assert actual.dtype == np.int64
	assert actual.ndim == 1
	assert actual.flags.c_contiguous


@pytest.mark.parametrize("invalid", [np.array([[0, 2]]), np.array([0.5]), np.array(0)])
def test_intervals_reject_invalid_degree_arrays(
	module: AugmentedPersistenceModule,
	invalid: np.ndarray,
) -> None:
	with pytest.raises(TypeError, match="degrees"):
		module.intervals(invalid)


def test_sparse_dataclass_equality_uses_array_values(sixpack_with_representatives: SixPack) -> None:
	pack = sixpack_with_representatives
	assert pack.representatives is not None
	assert pack.morphisms is not None
	for data in (pack.representatives, pack.morphisms):
		copied = type(data)(
			**{
				f.name: {b: row.copy() for b, row in getattr(data, f.name).items()}
				for f in fields(data)
			},
		)
		assert data == copied
		assert data != object()
		for f in fields(copied):
			rows = getattr(copied, f.name)
			if rows:
				next(iter(rows.values()))[0] += 1
				break
		assert data != copied
