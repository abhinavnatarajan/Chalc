"""Augmented persistence modules and 6-packs of persistence modules."""

from __future__ import annotations

from collections.abc import Callable, Iterator, Mapping, Sequence
from numbers import Integral
from typing import ClassVar, Literal, Self, get_args, overload

import h5py
import numpy as np

type PersistenceModuleName = Literal[
	"kernel",
	"cokernel",
	"domain",
	"codomain",
	"image",
	"relative",
]
"""Names of persistence modules in a 6-pack of persistence modules."""


class AugmentedPersistenceModule:
	r"""An :math:`\mathbb{R}`-indexed persistence module of :math:`\mathbb{N}`-graded vector spaces.

	The module is augmented with (1) a choice of interval basis, where each interval is a pairing
	between a generator and a relation, and (2) a strict linear ordering of the set of all
	generators and relations, which is compatible with their entrance times.

	In particular, each interval is uniquely identified by the integer index of its generator,
	and distinct integer indices may correspond to elements that have the same entrance time.
	"""

	__slots__ = ("_degrees", "_entrance_times", "_intervals")
	_intervals: dict[int, int | None]
	_degrees: dict[int, int]
	_entrance_times: np.ndarray[tuple[int], np.dtype[np.floating]]

	def __init__(
		self,
		intervals: Mapping[int, int | None] | None = None,
		entrance_times: Sequence[float] | np.ndarray[tuple[int], np.dtype[np.floating]] = (),
		degrees: Mapping[int, int] | None = None,
	) -> None:
		r"""Initialize the persistence module from its augmented interval decomposition.

		Args:
			intervals: A sequence of intervals, given as a mapping from generators to relations.
			entrance_times: The sorted sequence of entrance times of the generators and relations
				that define the persistence module.
			degrees: A sequence associating each interval (identified by its generator index)
				to its degree.
		"""
		self._intervals = dict(intervals) if intervals is not None else {}
		self._degrees = dict(degrees) if degrees is not None else {}
		self._entrance_times = np.array(entrance_times, dtype=np.float64, copy=True)
		if (
			self._entrance_times.ndim != 1
			or np.isnan(self._entrance_times).any()
			or np.any(self._entrance_times[:-1] > self._entrance_times[1:])
		):
			msg = "Endpoint entrance times must be a nondecreasing sequence without NaN."
			raise ValueError(msg)
		if self._degrees.keys() != self._intervals.keys():
			msg = "Degrees must be specified for exactly the interval generators."
			raise ValueError(msg)
		for birth, death in self._intervals.items():
			if (
				not isinstance(birth, Integral)
				or not 0 <= birth < len(self._entrance_times)
				or (
					death is not None
					and (
						not isinstance(death, Integral)
						or not birth <= death < len(self._entrance_times)
					)
				)
			):
				msg = "Interval indices must be in range, with birth preceding death."
				raise ValueError(msg)
			if not isinstance(self._degrees[birth], Integral) or self._degrees[birth] < 0:
				msg = "Interval degrees must be nonnegative integers."
				raise ValueError(msg)
		self._entrance_times.setflags(write=False)

	def __eq__(self, other: object) -> bool:
		"""Check if two persistence modules are equal."""
		if not isinstance(other, AugmentedPersistenceModule):
			return NotImplemented
		return (
			self._intervals == other._intervals
			and np.array_equal(self._entrance_times, other._entrance_times)
			and self._degrees == other._degrees
		)

	def __bool__(self) -> bool:
		"""Check if the persistence module is non-empty."""
		return bool(self._intervals)

	def __hash__(self) -> int:  # noqa: D105
		return hash(
			(
				tuple(self._entrance_times),
				frozenset(self._degrees.items()),
				frozenset(self._intervals.items()),
			),
		)

	@property
	def degrees(self) -> list[int]:
		"""Return all inhabited degrees in the persistence module."""
		return sorted(set(self._degrees.values()))

	@property
	def entrance_times(self) -> list[float]:
		"""Return all unique endpoint entrance times in the persistence module."""
		return np.unique(self._entrance_times).tolist()

	@overload
	def intervals(
		self,
		degrees: int,
	) -> np.ndarray[tuple[int, Literal[2]], np.dtype[np.float64]]: ...

	@overload
	def intervals(
		self,
		degrees: Sequence[int] | None = None,
	) -> dict[int, np.ndarray[tuple[int, Literal[2]], np.dtype[np.float64]]]: ...

	def intervals(
		self,
		degrees=None,
	):
		r"""
		Get the interval decomposition of the persistence module.

		Args:
			degrees:
				Degrees in which the intervals are desired.
				If a list is provided then a list of matrices is returned,
				with the order of matrices respecting the order of entries of `degree`.
				If `degree` is not provided then the returned matrix will contain
				intervals features from all degrees from zero to ``max(self.degrees)``.

		Returns:
			An :math:`m \times 2` matrix whose rows are the endpoints of intervals in the
			interval decomposition, or a dictionary of such matrices indexed by degree.
			The rows are sorted by generator index.
			Degenerate intervals (having zero persistence) are not removed.

		"""
		if degrees is None:
			degrees = self.degrees
		if isinstance(degrees, Integral):
			result = [
				[self._entrance_times[b], self._entrance_times[d] if d is not None else np.inf]
				for b, d in sorted(self._intervals.items())
				if self._degrees[b] == degrees
			]
			return np.asarray(result, dtype=np.float64).reshape((-1, 2))
		if isinstance(degrees, Sequence) and all(isinstance(d, Integral) for d in degrees):
			return {d: self.intervals(d) for d in degrees}
		msg = "`degrees` must be an integer or a sequence of integers."
		raise TypeError(msg)


class SixPack(Mapping[PersistenceModuleName, dict[int, int | None]]):
	"""6-pack of persistence modules."""

	__slots__ = (
		"_codomain",
		"_cokernel",
		"_degrees",
		"_domain",
		"_entrance_times",
		"_image",
		"_kernel",
		"_relative",
	)
	_degrees: np.ndarray[tuple[int], np.dtype[np.int64]]
	_entrance_times: np.ndarray[tuple[int], np.dtype[np.float64]]
	_kernel: dict[int, int | None]
	_cokernel: dict[int, int | None]
	_domain: dict[int, int | None]
	_codomain: dict[int, int | None]
	_image: dict[int, int | None]
	_relative: dict[int, int | None]
	__hash__ = None
	_MODULE_NAMES: ClassVar[tuple[PersistenceModuleName, ...]] = get_args(
		PersistenceModuleName.__value__,
	)

	def __init__(
		self,
		*,
		domain: Mapping[int, int | None] | None = None,
		codomain: Mapping[int, int | None] | None = None,
		kernel: Mapping[int, int | None] | None = None,
		image: Mapping[int, int | None] | None = None,
		cokernel: Mapping[int, int | None] | None = None,
		relative: Mapping[int, int | None] | None = None,
		entrance_times: Sequence[float] = (),
		degrees: Sequence[int] = (),
	) -> None:
		"""Initialise the 6-pack of persistence modules."""
		self._kernel = dict(kernel) if kernel is not None else {}
		self._cokernel = dict(cokernel) if cokernel is not None else {}
		self._domain = dict(domain) if domain is not None else {}
		self._codomain = dict(codomain) if codomain is not None else {}
		self._image = dict(image) if image is not None else {}
		self._relative = dict(relative) if relative is not None else {}
		self._entrance_times = np.array(entrance_times, dtype=np.float64, copy=True)
		self._degrees = np.array(degrees, dtype=np.int64, copy=True)
		if self._entrance_times.ndim != 1 or self._degrees.shape != self._entrance_times.shape:
			msg = "Entrance times and degrees must be one-dimensional arrays of equal length."
			raise ValueError(msg)

	@overload
	def get_module(
		self,
		key: PersistenceModuleName,
		*,
		return_index_map: Literal[False] = False,
	) -> AugmentedPersistenceModule:
		...

	@overload
	def get_module(
		self,
		key: PersistenceModuleName,
		*,
		return_index_map: Literal[True],
	) -> tuple[AugmentedPersistenceModule, dict[int, int]]:
		...

	@overload
	def get_module(
		self,
		key: PersistenceModuleName,
		*,
		return_index_map: bool,
	) -> AugmentedPersistenceModule | tuple[AugmentedPersistenceModule, dict[int, int]]:
		...

	def get_module(
		self,
		key: PersistenceModuleName,
		*,
		return_index_map: bool = False,
	) -> AugmentedPersistenceModule | tuple[AugmentedPersistenceModule, dict[int, int]]:
		"""Access a specific persistence module in the 6-pack.

		Constructs the persistence module from scratch on each access.
		Cache the results if performance is a concern.
		"""
		module = self[key]
		endpoint_indices = [*module.keys(), *(d for d in module.values() if d is not None)]
		sorted_endpoint_indices = sorted(endpoint_indices)
		new_indices = {e: i for i, e in enumerate(sorted_endpoint_indices)}
		intervals = {
			new_indices[b]: new_indices[d] if d is not None else None for b, d in module.items()
		}
		entrance_times = self._entrance_times[sorted_endpoint_indices]
		degrees = {new_indices[b]: int(self._degrees[b]) - (key == "kernel") for b in module}
		if return_index_map:
			return AugmentedPersistenceModule(intervals, entrance_times, degrees), new_indices
		return AugmentedPersistenceModule(intervals, entrance_times, degrees)

	def __getitem__(self, key: object) -> dict[int, int | None]:
		"""Get a persistence module as raw birth-death integer indices.

		In most cases, prefer :func:`SixPack.get_module`.
		"""
		if isinstance(key, str) and key in self._MODULE_NAMES:
			return getattr(self, "_" + key)
		raise KeyError(key)

	def __setitem__(self, key: PersistenceModuleName, value: Mapping[int, int | None]) -> None:
		"""Set pairings for a persistence module in the 6-pack, for internal use only."""
		if key not in self._MODULE_NAMES:
			raise KeyError(key)
		setattr(self, "_" + key, dict(value))

	def __iter__(self) -> Iterator[PersistenceModuleName]:
		"""Iterate over the six persistence module names."""
		return iter(self._MODULE_NAMES)

	def __len__(self) -> int:
		"""Get the number of persistence modules in the 6-pack (hint: it's in the name).

		Implemented for conformance with the Mapping protocol.
		"""
		return len(self._MODULE_NAMES)

	def __bool__(self) -> bool:
		"""Check whether any persistence module in the 6-pack is non-trivial."""
		return any(self.values())

	def __eq__(self, other: object) -> bool:
		"""Check equality of 6-packs.

		Checks if all the persistence modules in the 6-packs are equal,
		with identical interval endpoint metadata.
		"""
		if not isinstance(other, SixPack):
			return NotImplemented
		return (
			all(self[name] == other[name] for name in self)
			and np.array_equal(self._entrance_times, other._entrance_times)
			and np.array_equal(self._degrees, other._degrees)
		)

	@classmethod
	def module_names(cls) -> tuple[PersistenceModuleName, ...]:
		"""Return all six persistence module names in iteration order."""
		return cls._MODULE_NAMES

	def max_nonempty_degree(self) -> int:
		"""Get the maximum degree of all features in all persistence modules in the 6-pack.

		Returns -1 if all persistence modules in the 6-pack are trivial.
		"""
		return max(
			(
				int(self._degrees[birth_idx]) - (name == "kernel")
				for name in self
				for birth_idx in self[name]
			),
			default=-1,
		)

	def filter(
		self,
		predicate: Callable[[PersistenceModuleName, int, int, int | None, float, float], bool],
	) -> Self:
		"""Retain intervals satisfying ``predicate``, modifying this 6-pack in place.

		The predicate receives the persistence module name, degree, birth index, death index,
		birth time, and death time of a feature. Infinite intervals have death index
		``None`` and death time ``np.inf``.

		Shared endpoint metadata is retained. Returns this 6-pack for chaining.
		"""
		for name in self:
			self[name] = {
				b: d
				for b, d in self[name].items()
				if predicate(
					name,
					int(self._degrees[b]) - (name == "kernel"),
					b,
					d,
					float(self._entrance_times[b]),
					float(self._entrance_times[d]) if d is not None else np.inf,
				)
			}
		return self

	def cleanup(self) -> Self:
		"""Remove unused endpoints and reindex all persistence modules in the 6-pack.

		This is an in-place operation. Entrance times, degrees, and intervals are preserved.
		This is useful after :meth:`filter`.

		Returns this 6-pack for method chaining.
		"""
		active_indices = {birth_idx for module in self.values() for birth_idx in module}
		active_indices.update(
			death_idx
			for module in self.values()
			for death_idx in module.values()
			if death_idx is not None
		)
		sorted_indices = sorted(active_indices)
		index_map = {old: new for new, old in enumerate(sorted_indices)}
		self._degrees = self._degrees[sorted_indices]
		self._entrance_times = self._entrance_times[sorted_indices]
		for name in self:
			self[name] = {
				index_map[b]: index_map[d] if d is not None else None for b, d in self[name].items()
			}
		return self

	def to_hdf5(self, file: h5py.Group) -> None:
		"""
		Save a 6-pack of persistence modules to a HDF5 file or group.

		Args:
			file : A h5py file or group.

		"""
		file["entrance_times"] = self._entrance_times
		file["degrees"] = self._degrees
		for name in self:
			file[name] = np.asarray(
				[(b, d if d is not None else -1) for b, d in self[name].items()],
				dtype=np.int64,
			).reshape((-1, 2))

	@classmethod
	def from_hdf5(cls, file: h5py.Group) -> Self:
		"""Load a 6-pack written by :meth:`to_hdf5` from an HDF5 file or group.

		Raises ``RuntimeError`` if a required dataset is missing or has an invalid
		shape or dtype.
		"""
		fields = ["entrance_times", "degrees", *cls._MODULE_NAMES]
		missing = [field for field in fields if field not in file]
		if missing:
			msg = f"Invalid file: missing fields {missing}"
			raise RuntimeError(msg)
		arrays = {}
		for field in fields:
			dataset = file[field]
			shape_valid = (
				(
					dataset.ndim == 1
					if field in {"entrance_times", "degrees"}
					else dataset.ndim == 2 and dataset.shape[1] == 2  # noqa: PLR2004
				)
				if isinstance(dataset, h5py.Dataset)
				else False
			)
			dtype = np.float64 if field == "entrance_times" else np.int64
			if not shape_valid or dataset.dtype != dtype:
				msg = f"Invalid dataset: {field}"
				raise RuntimeError(msg)
			arrays[field] = dataset[:]
		if arrays["degrees"].shape != arrays["entrance_times"].shape:
			msg = "Invalid file: entrance times and degrees have different lengths."
			raise RuntimeError(msg)
		modules = {
			name: {int(b): int(d) if d != -1 else None for b, d in arrays[name]}
			for name in cls._MODULE_NAMES
		}
		return cls(**modules, entrance_times=arrays["entrance_times"], degrees=arrays["degrees"])
