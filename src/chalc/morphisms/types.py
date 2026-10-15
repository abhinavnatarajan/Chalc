"""Augmented persistence modules and 6-packs of persistence modules."""

from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, field, fields
from numbers import Integral
from types import MappingProxyType
from typing import ClassVar, Literal, get_args, overload

import h5py
import numpy as np

__all__ = [
	"AugmentedPersistenceModule",
	"PersistenceModuleName",
	"SixPack",
	"SixPackMorphisms",
	"SixPackPairings",
	"SixPackRepresentatives",
]

type PersistenceModuleName = Literal[
	"kernel",
	"domain",
	"image",
	"codomain",
	"cokernel",
	"relative",
]
"""Names of persistence modules in a 6-pack of persistence modules."""


class AugmentedPersistenceModule:
	r"""An :math:`\mathbb{R}`-indexed persistence module of :math:`\mathbb{N}`-graded vector spaces.

	The module is augmented with (1) a choice of interval basis, where each interval is a pairing
	between a generator and a relation, (2) a strict linear ordering of the set of all
	generators and relations, which is compatible with their entrance times.

	In particular, each interval is uniquely identified by the integer index of its generator,
	and distinct integer indices may correspond to elements that have the same entrance time.
	"""

	__slots__ = ("_degrees", "_entrance_times", "_intervals")
	_intervals: dict[int, int | None]
	_degrees: dict[int, int]
	_entrance_times: np.ndarray[tuple[int], np.dtype[np.float64]]

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

	def __len__(self) -> int:
		"""Return the number of features in the module."""
		return len(self._intervals)

	def nonempty_degrees(self) -> list[int]:
		"""Return all inhabited degrees in the persistence module."""
		return np.unique(list(self._degrees.values())).tolist()

	@property
	def pairings(self) -> MappingProxyType[int, int | None]:
		"""Get the birth-death pairing indices for this module."""
		return MappingProxyType(self._intervals)

	@property
	def entrance_times(self) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
		"""Return all entrance times of generators/relations in the persistence module.

		The entrance times are sorted by generator/relation index.
		"""
		return self._entrance_times

	def interval(self, i: int) -> tuple[float, float]:
		"""Return the interval with birth index ``i``."""
		if i not in self._intervals:
			raise KeyError
		death_idx = self._intervals[i]
		return (
			self._entrance_times[i],
			self._entrance_times[death_idx] if death_idx is not None else float("inf"),
		)

	@overload
	def intervals(
		self,
		degrees: int,
	) -> np.ndarray[tuple[int, Literal[2]], np.dtype[np.float64]]: ...

	@overload
	def intervals(
		self,
		degrees: Sequence[int] | None = ...,
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
			degrees = self.nonempty_degrees()
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

	def to_hdf5(
		self,
		file: h5py.Group,
		*,
		pairings: Mapping[int, int | None] | None = None,
	) -> None:
		"""Save interval pairings, endpoint entrance times, and interval degrees.

		Infinite deaths are stored as ``-1``. Degrees are stored in pairing row order.
		By default, pairings use this module's endpoint indices. An enclosing 6-pack
		can supply its original pairings instead, in the same interval order.
		"""
		stored_pairings = self._intervals if pairings is None else pairings
		file["pairings"] = np.asarray(
			[
				(birth, death if death is not None else -1)
				for birth, death in stored_pairings.items()
			],
			dtype=np.int64,
		).reshape((-1, 2))
		file["entrance_times"] = self._entrance_times
		file["degrees"] = np.asarray(
			[self._degrees[birth] for birth in self._intervals],
			dtype=np.int64,
		)


class SixPack:
	"""6-pack of persistence modules."""

	__slots__ = (
		"_codomain",
		"_cokernel",
		"_domain",
		"_endpoint_relabelling",
		"_image",
		"_kernel",
		"_morphisms",
		"_pairings",
		"_relative",
		"_representatives",
		"metadata",
	)
	_codomain: AugmentedPersistenceModule
	_cokernel: AugmentedPersistenceModule
	_domain: AugmentedPersistenceModule
	_image: AugmentedPersistenceModule
	_kernel: AugmentedPersistenceModule
	_relative: AugmentedPersistenceModule
	# Pairings in original ambient coordinates.
	_pairings: SixPackPairings
	# Original endpoint indices to relabelled coordinates, after filtering
	# out degenerate intervals and removing gaps in the numbering, per module.
	_endpoint_relabelling: dict[PersistenceModuleName, dict[int, int]]
	_morphisms: SixPackMorphisms | None
	_representatives: SixPackRepresentatives | None
	metadata: dict | None  # Optional metadata
	__hash__ = None
	_MODULE_NAMES: ClassVar[tuple[PersistenceModuleName, ...]] = get_args(
		PersistenceModuleName.__value__,
	)

	def __init__(
		self,
		*,
		pairings: SixPackPairings | None = None,
		entrance_times: Sequence[float] | np.ndarray[tuple[int], np.dtype[np.floating]] = (),
		degrees: Sequence[int] | np.ndarray[tuple[int], np.dtype[np.integer]] = (),
		representatives: SixPackRepresentatives | None = None,
		morphisms: SixPackMorphisms | None = None,
	) -> None:
		"""Initialise the 6-pack of persistence modules.

		Assuming a linear ordering of all interval endpoints across
		all diagrams, that is compatible with their entrance times,
		each module should be provided as a collection of pairings
		between endpoint indices (or ``birth index -> None`` for
		infinite intervals). Entrance times and degrees should
		correspond to the chosen ordering. The supplied frozen dataclasses are
		owned by this 6-pack; representatives and maps are reindexed in place.
		Consumer metadata can be assigned to ``metadata`` after construction.
		"""
		self._pairings = pairings if pairings is not None else SixPackPairings()
		self._representatives = representatives
		self._morphisms = morphisms
		self._endpoint_relabelling = {}
		self.metadata = None

		et = np.asarray(entrance_times, dtype=np.float64)
		degs = np.asarray(degrees, dtype=np.int64)
		if et.ndim != 1 or degs.shape != et.shape:
			msg = "Entrance times and degrees must be one-dimensional arrays of equal length."
			raise ValueError(msg)
		if np.isnan(et).any() or np.any(et[:-1] > et[1:]):
			msg = "Endpoint entrance times must be a nondecreasing sequence without NaN."
			raise ValueError(msg)
		self._init_modules_and_reps(et, degs)

	def _init_modules_and_reps(
		self,
		entrance_times: np.ndarray[tuple[int], np.dtype[np.floating]],
		degrees: Sequence[int] | np.ndarray[tuple[int], np.dtype[np.integer]],
	) -> None:
		# Keep the original generator sets for validation before filtering the
		# owned pairing dictionaries and reindexing representatives and maps.
		original_pairings = {name: getattr(self._pairings, name) for name in self._MODULE_NAMES}
		for name in self._MODULE_NAMES:
			module = self._init_module(name, entrance_times, degrees)
			setattr(self, "_" + name, module)
		if self._representatives is not None:
			for name in self._MODULE_NAMES:
				reps = getattr(self._representatives, name)
				if not reps.keys() <= original_pairings[name].keys():
					msg = f"Representatives for {name} must refer to interval generators."
					raise ValueError(msg)
				idx_map = self._endpoint_relabelling[name]
				object.__setattr__(
					self._representatives,
					name,
					{
						idx_map[birth]: support
						for birth, support in reps.items()
						if birth in idx_map
					},
				)
		if self._morphisms is not None:
			for morph_field in fields(self._morphisms):
				name = morph_field.name
				source, target = _morphism_modules(name)
				morphism = getattr(self._morphisms, name)
				if not morphism.keys() <= original_pairings[source].keys() or any(
					not set(support) <= original_pairings[target].keys()
					for support in morphism.values()
				):
					msg = f"Map {name} must refer to source and target interval generators."
					raise ValueError(msg)
				source_indices = self._endpoint_relabelling[source]
				target_indices = self._endpoint_relabelling[target]
				object.__setattr__(
					self._morphisms,
					name,
					{
						source_indices[birth]: [
							target_indices[idx] for idx in support if idx in target_indices
						]
						for birth, support in morphism.items()
						if birth in source_indices
					},
				)

	def _init_module(
		self,
		key: PersistenceModuleName,
		entrance_times: np.ndarray[tuple[int], np.dtype[np.floating]],
		degrees: Sequence[int] | np.ndarray[tuple[int], np.dtype[np.integer]],
	) -> AugmentedPersistenceModule:
		"""Initialize persistence modules from pairing data.

		This function reindexes the interval endpoints so that they are numbered sequentially
		without gaps; the numbering of endpoints may differ from the raw birth-death indices
		since not all endpoints in the 6-pack may appear in the module. Degenerate intervals
		are discarded.
		"""
		module_pairings = getattr(self._pairings, key)
		for birth, death in module_pairings.items():
			if (
				not isinstance(birth, Integral)
				or not 0 <= int(birth) < len(entrance_times)
				or (
					death is not None
					and (
						not isinstance(death, Integral)
						or not int(birth) <= int(death) < len(entrance_times)
					)
				)
			):
				msg = f"Invalid pairings for {key}: indices must be in range and ordered."
				raise ValueError(msg)
			if int(degrees[int(birth)]) - (key == "kernel") < 0:
				msg = "Interval degrees must be nonnegative integers."
				raise ValueError(msg)
		# Discard degenerate intervals, retaining the original indices in pairings.
		filtered_module_pairings = {
			b: d
			for b, d in module_pairings.items()
			if entrance_times[b] < (np.inf if d is None else entrance_times[d])
		}
		object.__setattr__(self._pairings, key, filtered_module_pairings)
		sorted_endpoint_indices = _endpoint_indices(filtered_module_pairings)
		new_indices = {e: i for i, e in enumerate(sorted_endpoint_indices)}
		intervals = {
			new_indices[b]: new_indices[d] if d is not None else None
			for b, d in filtered_module_pairings.items()
		}
		endpoint_entrance_times = entrance_times[sorted_endpoint_indices]
		interval_degrees = {
			new_indices[interval_birth_idx]: int(degrees[interval_birth_idx]) - (key == "kernel")
			for interval_birth_idx in filtered_module_pairings
		}
		self._endpoint_relabelling[key] = {
			birth: new_indices[birth] for birth in filtered_module_pairings
		}
		return AugmentedPersistenceModule(
			intervals,
			endpoint_entrance_times,
			interval_degrees,
		)

	def __bool__(self) -> bool:
		"""Check whether any persistence module in the 6-pack is non-trivial."""
		return any(module for _, module in self)

	def __eq__(self, other: object) -> bool:
		"""Check equality of 6-packs.

		Checks if all the persistence modules in the 6-packs are equal,
		with identical interval endpoint metadata.
		"""
		if not isinstance(other, SixPack):
			return NotImplemented
		return all(getattr(self, name) == getattr(other, name) for name, _ in self)

	@property
	def kernel(self) -> AugmentedPersistenceModule:
		"""Kernel persistence module."""
		return self._kernel

	@property
	def domain(self) -> AugmentedPersistenceModule:
		"""Domain persistence module."""
		return self._domain

	@property
	def image(self) -> AugmentedPersistenceModule:
		"""Image persistence module."""
		return self._image

	@property
	def codomain(self) -> AugmentedPersistenceModule:
		"""Codomain persistence module."""
		return self._codomain

	@property
	def cokernel(self) -> AugmentedPersistenceModule:
		"""Cokernel persistence module."""
		return self._cokernel

	@property
	def relative(self) -> AugmentedPersistenceModule:
		"""Relative persistence module."""
		return self._relative

	@property
	def pairings(self) -> SixPackPairings:
		"""Birth-death pairings defining the 6-pack."""
		return self._pairings

	@property
	def morphisms(self) -> SixPackMorphisms | None:
		"""Natural morphisms between modules in the 6-pack."""
		return self._morphisms

	@property
	def representatives(self) -> SixPackRepresentatives | None:
		"""Representative cycles for intervals in the modules in the 6-pack."""
		return self._representatives

	@classmethod
	def module_names(cls) -> tuple[PersistenceModuleName, ...]:
		"""Return all six persistence module names in iteration order."""
		return cls._MODULE_NAMES

	def __iter__(self) -> Iterator[tuple[PersistenceModuleName, AugmentedPersistenceModule]]:
		"""Iterate over the persistence modules in the 6-pack."""
		return ((name, getattr(self, name)) for name in self._MODULE_NAMES)

	def max_nonempty_degree(self) -> int:
		"""Get the maximum degree of all features in all persistence modules in the 6-pack.

		Returns -1 if all persistence modules in the 6-pack are trivial.
		"""
		return max(
			(degree for _, module in self for degree in module.nonempty_degrees()),
			default=-1,
		)

	def to_hdf5(self, file: h5py.Group) -> None:
		"""Save modules, raw pairings, and optional representatives and maps.

		Each module stores its own endpoint times and interval degrees. Representatives
		and maps use the same module-local indices as in memory. Consumer metadata is
		not serialized; consumers may store it separately in the same file or group.
		"""
		for name, module in self:
			module.to_hdf5(file.create_group(name), pairings=getattr(self._pairings, name))
		if self._representatives is not None:
			group = file.create_group("representatives")
			for name in self._MODULE_NAMES:
				_write_sparse_rows(group.create_group(name), getattr(self._representatives, name))
		if self._morphisms is not None:
			group = file.create_group("morphisms")
			for morph_field in fields(self._morphisms):
				_write_sparse_rows(
					group.create_group(morph_field.name),
					getattr(self._morphisms, morph_field.name),
				)

	@classmethod
	def from_hdf5(cls, file: h5py.Group) -> SixPack:
		"""Load modules, raw pairings, and optional representatives and maps.

		Consumer metadata is ignored and initialized to ``None``.
		Malformed datasets raise ``RuntimeError``.
		"""
		try:
			representatives = None
			if "representatives" in file:
				group = _read_group(file, "representatives")
				representatives = SixPackRepresentatives(
					**{
						name: _read_sparse_rows(_read_group(group, name))
						for name in cls._MODULE_NAMES
					},
				)
			morphisms = None
			if "morphisms" in file:
				group = _read_group(file, "morphisms")
				morphisms = SixPackMorphisms(
					**{
						f.name: _read_sparse_rows(_read_group(group, f.name))
						for f in fields(SixPackMorphisms)
					},
				)
			self = cls.__new__(cls)
			self._endpoint_relabelling = {}
			pairings_dict = {}
			for name in cls._MODULE_NAMES:
				group = _read_group(file, name)
				pairings = _read_pairings(group, "pairings")
				module, self._endpoint_relabelling[name] = _read_module(group, pairings)
				setattr(self, "_" + name, module)
				pairings_dict[name] = pairings
			self._pairings = SixPackPairings(**pairings_dict)
			self._representatives = representatives
			self._morphisms = morphisms
			self.metadata = None
			self._validate_local_generators()
		except (KeyError, IndexError, TypeError, ValueError, OverflowError) as err:
			msg = f"Invalid HDF5 SixPack: {err}"
			raise RuntimeError(msg) from err
		else:
			return self

	def _validate_local_generators(self) -> None:
		if self._representatives is not None:
			for name, module in self:
				if not getattr(self._representatives, name).keys() <= module.pairings.keys():
					msg = f"Representatives for {name} must refer to interval generators."
					raise ValueError(msg)
		if self._morphisms is not None:
			for morph_field in fields(self._morphisms):
				name = morph_field.name
				source, target = _morphism_modules(name)
				mapping = getattr(self._morphisms, name)
				if not mapping.keys() <= getattr(self, source).pairings.keys() or any(
					not set(support) <= getattr(self, target).pairings.keys()
					for support in mapping.values()
				):
					msg = f"Map {name} must refer to source and target interval generators."
					raise ValueError(msg)


@dataclass(slots=True, frozen=True)
class SixPackMorphisms:
	"""Natural maps between persistence modules in a 6-pack."""

	kernel_to_domain: dict[int, list[int]] = field(kw_only=True)
	"""Monomorphism of kernel into the domain."""

	domain_to_image: dict[int, list[int]] = field(kw_only=True)
	"""Epimorphism of domain onto image."""

	image_to_codomain: dict[int, list[int]] = field(kw_only=True)
	"""Monomorphism of the image into the codomain."""

	domain_to_codomain: dict[int, list[int]] = field(kw_only=True)
	"""Domain to codomain map."""

	codomain_to_cokernel: dict[int, list[int]] = field(kw_only=True)
	"""Epimorphism of codomain onto cokernel."""


@dataclass(slots=True, frozen=True)
class SixPackRepresentatives:
	"""Representative cycles for persistence modules in a 6-pack."""

	kernel: dict[int, list[int]] = field(kw_only=True)
	"""Kernel representatives."""

	domain: dict[int, list[int]] = field(kw_only=True)
	"""Domain representatives."""

	image: dict[int, list[int]] = field(kw_only=True)
	"""Image representatives."""

	codomain: dict[int, list[int]] = field(kw_only=True)
	"""Codomain representatives."""

	cokernel: dict[int, list[int]] = field(kw_only=True)
	"""Cokernel representatives."""

	relative: dict[int, list[int]] = field(kw_only=True)
	"""Relative representatives."""


@dataclass(slots=True, frozen=True)
class SixPackPairings:
	"""Birth-death pairings for persistence modules in a 6-pack."""

	kernel: dict[int, int | None] = field(default_factory=dict, kw_only=True)
	"""Kernel pairings."""

	domain: dict[int, int | None] = field(default_factory=dict, kw_only=True)
	"""Domain pairings."""

	image: dict[int, int | None] = field(default_factory=dict, kw_only=True)
	"""Image pairings."""

	codomain: dict[int, int | None] = field(default_factory=dict, kw_only=True)
	"""Codomain pairings."""

	cokernel: dict[int, int | None] = field(default_factory=dict, kw_only=True)
	"""Cokernel pairings."""

	relative: dict[int, int | None] = field(default_factory=dict, kw_only=True)
	"""Relative pairings."""


def _endpoint_indices(pairings: Mapping[int, int | None]) -> list[int]:
	return sorted(set(pairings) | {death for death in pairings.values() if death is not None})


def _morphism_modules(name: str) -> tuple[PersistenceModuleName, PersistenceModuleName]:
	source, target = name.split("_to_")
	module_names: dict[str, PersistenceModuleName] = {name: name for name in SixPack.module_names()}
	return module_names[source], module_names[target]


def _read_pairings(group: h5py.Group, name: str) -> dict[int, int | None]:
	rows = _read_array(group, name, np.dtype(np.int64), 2)
	if (
		rows.shape[1] != 2  # noqa: PLR2004
		or len(np.unique(rows[:, 0])) != len(rows)
		or np.any(rows[:, 0] < 0)
		or np.any((rows[:, 1] != -1) & (rows[:, 1] < rows[:, 0]))
	):
		msg = f"Invalid pairings dataset: {name}"
		raise ValueError(msg)
	return {int(birth): int(death) if death != -1 else None for birth, death in rows}


def _read_module(
	group: h5py.Group,
	pairings: Mapping[int, int | None],
) -> tuple[AugmentedPersistenceModule, dict[int, int]]:
	indices = _endpoint_indices(pairings)
	times = _read_array(group, "entrance_times", np.dtype(np.float64), 1)
	degrees = _read_array(group, "degrees", np.dtype(np.int64), 1)
	if len(times) != len(indices) or len(degrees) != len(pairings):
		msg = f"Invalid endpoint or degree count: {group.name}"
		raise ValueError(msg)
	index_map = {old: new for new, old in enumerate(indices)}
	module = AugmentedPersistenceModule(
		{index_map[b]: index_map[d] if d is not None else None for b, d in pairings.items()},
		times,
		{index_map[b]: int(degree) for b, degree in zip(pairings, degrees, strict=True)},
	)
	if any(
		times[birth] == (np.inf if death is None else times[death])
		for birth, death in module.pairings.items()
	):
		msg = f"Degenerate intervals are not allowed in a SixPack: {group.name}"
		raise ValueError(msg)
	return module, {birth: index_map[birth] for birth in pairings}


def _write_sparse_rows(group: h5py.Group, rows: Mapping[int, Sequence[int]]) -> None:
	group["keys"] = np.asarray(list(rows), dtype=np.int64)
	group["offsets"] = np.concatenate(
		(
			np.zeros(1, dtype=np.int64),
			np.cumsum([len(row) for row in rows.values()], dtype=np.int64),
		),
	)
	group["indices"] = np.asarray([idx for row in rows.values() for idx in row], dtype=np.int64)


def _read_sparse_rows(group: h5py.Group) -> dict[int, list[int]]:
	keys = _read_array(group, "keys", np.dtype(np.int64), 1)
	offsets = _read_array(group, "offsets", np.dtype(np.int64), 1)
	indices = _read_array(group, "indices", np.dtype(np.int64), 1)
	if (
		len(offsets) != len(keys) + 1
		or offsets[0] != 0
		or offsets[-1] != len(indices)
		or np.any(offsets[:-1] > offsets[1:])
		or np.any(keys < 0)
		or np.any(indices < 0)
		or len(np.unique(keys)) != len(keys)
	):
		msg = f"Invalid sparse rows: {group.name}"
		raise ValueError(msg)
	return {int(key): indices[offsets[i] : offsets[i + 1]].tolist() for i, key in enumerate(keys)}


def _read_group(group: h5py.Group, name: str) -> h5py.Group:
	child = group[name]
	if not isinstance(child, h5py.Group):
		msg = f"Invalid group: {name}"
		raise TypeError(msg)
	return child


def _read_array(group: h5py.Group, name: str, dtype: np.dtype, ndim: int) -> np.ndarray:
	dataset = group[name]
	if not isinstance(dataset, h5py.Dataset) or dataset.dtype != dtype or dataset.ndim != ndim:
		msg = f"Invalid dataset: {name}"
		raise ValueError(msg)
	return dataset[:]
