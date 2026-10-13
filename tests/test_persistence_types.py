"""Regression tests for augmented persistence modules and 6-packs."""

from pathlib import Path

import h5py
import numpy as np
import pytest

from chalc.morphisms import AugmentedPersistenceModule, SixPack


@pytest.fixture
def module() -> AugmentedPersistenceModule:
	return AugmentedPersistenceModule(
		{0: 2, 1: None, 3: 4},
		[0.0, 0.0, 1.0, 2.0, 3.0],
		{0: 0, 1: 2, 3: 2},
	)


def test_intervals_by_degree(module: AugmentedPersistenceModule) -> None:
	assert module.degrees == [0, 2]
	assert module.entrance_times == [0.0, 1.0, 2.0, 3.0]
	np.testing.assert_array_equal(module.intervals(0), [[0.0, 1.0]])
	np.testing.assert_array_equal(module.intervals(2), [[0.0, np.inf], [2.0, 3.0]])
	assert module.intervals(1).shape == (0, 2)
	assert module.intervals(1).dtype == np.float64
	assert list(module.intervals()) == [0, 2]
	assert list(module.intervals([2, 1, 0])) == [2, 1, 0]
	# pyrefly: ignore [no-matching-overload]
	np.testing.assert_array_equal(module.intervals(np.int64(2)), module.intervals(2))


def test_empty_module() -> None:
	module = AugmentedPersistenceModule()
	assert not module
	assert module.degrees == []
	assert module.entrance_times == []
	assert module.intervals() == {}
	assert module.intervals(0).shape == (0, 2)


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


@pytest.fixture
def sixpack() -> SixPack:
	return SixPack(
		domain={1: 4, 2: None},
		kernel={3: 5},
		relative={3: None},
		entrance_times=[0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0],
		degrees=[0, 0, 0, 1, 1, 2, 3],
	)


def test_sixpack_mapping_and_modules(sixpack: SixPack) -> None:
	assert set(SixPack.module_names()) == {
		"kernel",
		"cokernel",
		"domain",
		"codomain",
		"image",
		"relative",
	}
	assert len(sixpack) == 6
	assert list(sixpack) == list(SixPack.module_names())
	assert len(list(sixpack.items())) == 6
	assert sixpack
	assert not SixPack()
	assert SixPack().max_nonempty_degree() == -1
	assert sixpack.max_nonempty_degree() == 1
	assert sixpack["domain"] == {1: 4, 2: None}
	assert sixpack.get_module("domain").degrees == [0]
	np.testing.assert_array_equal(
		sixpack.get_module("domain").intervals(0),
		[[0.0, 2.0], [0.0, np.inf]],
	)
	assert sixpack.get_module("kernel").degrees == [0]
	np.testing.assert_array_equal(sixpack.get_module("kernel").intervals(0), [[1.0, 3.0]])
	assert not sixpack.get_module("image")
	for key in ["ker", "unknown", None]:
		with pytest.raises(KeyError):
			sixpack[key]
		with pytest.raises(KeyError):
			# pyrefly: ignore [no-matching-overload]
			sixpack.get_module(key)
	with pytest.raises(TypeError):
		hash(sixpack)


def test_filter_and_cleanup_preserve_modules(sixpack: SixPack) -> None:
	calls = []

	def retain(  # noqa: PLR0917
		name: str,
		degree: int,
		b: int,
		d: int | None,
		birth: float,
		death: float,
	) -> bool:
		calls.append((name, degree, b, d, birth, death))
		return name != "relative"

	assert sixpack.filter(retain) is sixpack
	assert ("kernel", 0, 3, 5, 1.0, 3.0) in calls
	assert ("domain", 0, 2, None, 0.0, np.inf) in calls
	modules = {name: sixpack.get_module(name) for name in sixpack}
	assert sixpack.cleanup() is sixpack
	assert sixpack["domain"] == {0: 3, 1: None}
	assert sixpack["kernel"] == {2: 4}
	assert {name: sixpack.get_module(name) for name in sixpack} == modules
	assert sixpack.cleanup() is sixpack
	assert {name: sixpack.get_module(name) for name in sixpack} == modules
	assert sixpack.filter(lambda *_: False).cleanup() == SixPack()


def test_sixpack_equality(sixpack: SixPack) -> None:
	assert sixpack != object()
	assert SixPack(entrance_times=[0.0], degrees=[0]) != SixPack()
	assert SixPack(entrance_times=[0.0], degrees=[0]) != SixPack(entrance_times=[1.0], degrees=[0])
	assert SixPack(entrance_times=[0.0], degrees=[0]) != SixPack(entrance_times=[0.0], degrees=[1])


@pytest.mark.parametrize("empty", [False, True])
def test_hdf5_roundtrip(tmp_path: Path, sixpack: SixPack, *, empty: bool) -> None:
	original = SixPack() if empty else sixpack
	with h5py.File(tmp_path / "sixpack.h5", "w") as file:
		group = file.create_group("sixpack")
		original.to_hdf5(group)
		assert group["image"].shape == (0, 2)
		assert group["image"].dtype == np.int64
	with h5py.File(tmp_path / "sixpack.h5", "r") as file:
		loaded = SixPack.from_hdf5(file["sixpack"])
	assert loaded == original
	assert {name: loaded.get_module(name) for name in loaded} == {
		name: original.get_module(name) for name in original
	}


@pytest.mark.parametrize("damage", ["missing", "group", "dtype", "shape"])
def test_hdf5_invalid_dataset(tmp_path: Path, damage: str) -> None:
	with h5py.File(tmp_path / "invalid.h5", "w") as file:
		SixPack().to_hdf5(file)
		if damage == "missing":
			del file["degrees"]
		elif damage == "group":
			del file["domain"]
			file.create_group("domain")
		elif damage == "dtype":
			del file["domain"]
			file["domain"] = np.empty((0, 2), dtype=np.float64)
		else:
			del file["domain"]
			file["domain"] = np.empty((0,), dtype=np.int64)
		with pytest.raises(RuntimeError, match="Invalid"):
			SixPack.from_hdf5(file)
