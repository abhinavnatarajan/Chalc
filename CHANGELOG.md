# Changelog

## 15.0.0 (unreleased)

Changes since 14.2.0.

### Breaking API changes

- Move persistence types and filtration morphisms from `chalc.sixpack` to
  `chalc.morphisms`.
- Replace `SimplexPairings` with `AugmentedPersistenceModule`, which associates
  interval pairings with endpoint entrance times and homological degrees. Use
  `module.intervals(degrees)` to obtain birth/death arrays: an integer selects
  one array, while a sequence or `None` returns a dictionary indexed by degree.
  `module.pairings` is a read-only mapping from birth indices to death indices,
  with `None` representing infinite death.
- Access the modules of a `SixPack` through `kernel`, `domain`, `image`,
  `codomain`, `cokernel`, and `relative` attributes. The abbreviated diagram
  names (`ker`, `dom`, `im`, `cod`, `cok`, `rel`) are replaced by these full names
  in `PersistenceModuleName`. Iterating over a pack yields `(name, module)` pairs
  instead of mapping keys; the mapping interface and in-place `filter` API are
  removed.
- Construct a `SixPack` with `SixPackPairings`, entrance times, and simplex
  degrees, plus optional `SixPackRepresentatives` and `SixPackMorphisms`.
  Each module has compact endpoint indices; `pack.pairings` retains the original
  ambient indices. Distinct endpoints with equal entrance times remain distinct.
  The pack takes ownership of the supplied dataclasses and reindexes their data.
- Remove zero-persistence intervals when constructing a `SixPack`, including
  their representatives and coefficients in the natural maps.
- Rename `sixpack(max_diagram_dimension=...)` to
  `sixpack(max_homological_degree=...)`, `SixPack.names()` to
  `SixPack.module_names()`, and `max_nonempty_dimension()` to
  `max_nonempty_degree()`. Replace `pack.get_matrix(name, degree)` with
  `getattr(pack, name).intervals(degree)`. Entrance times and degrees now belong
  to individual modules rather than being exposed as shared arrays on the pack.
- Replace `SixPack.save` / `from_file` with `to_hdf5` / `from_hdf5` and a new HDF5
  layout. Files written by 14.2.0 must be converted or recomputed before loading
  with the new API.
- Replace `plot_diagram(pack, name, ...)` with
  `plot_persistence_diagram(module, ...)`. Plotting functions use `degrees`
  instead of `dimensions`. Optional arguments to `draw_filtration`, and
  `animation_length` and subsequent arguments to `animate_filtration`, are now
  keyword-only.
- Require `phimaker >= 0.4.0`.

### Added

- Add `sixpack(get_representatives=True)` to compute representative cycles for
  all six modules and the natural maps between kernel, domain, image, codomain,
  and cokernel. Access them through `pack.representatives` and `pack.morphisms`.
  Generator indices are local to each module; representative supports use
  ambient boundary-matrix indices, or mapping-cylinder indices for quotients.
- Preserve modules, original pairings, representatives, maps, and the endpoint
  relabelling through HDF5 round trips. Add `AugmentedPersistenceModule.to_hdf5`
  for writing individual modules. Consumer `pack.metadata` is assigned directly
  and is not serialized; quotient computations populate it with the maps from
  domain/codomain simplex indices into the mapping cylinder.
- Expose separate feature and truncation legend controls in
  `plot_persistence_diagram`, which consistently returns its axes, including for
  empty diagrams. `plot_sixpack` calculates default truncations for module names
  omitted from a supplied mapping.

### Fixed

- Apply homological degree limits correctly to kernel intervals and quotient
  relative homology, accounting for shifted degrees in the mapping cylinder.
- Validate interval indices, entrance-time ordering, and degrees before
  filtering; reject representatives and maps that refer to invalid generators.
  Reject malformed HDF5 datasets with a `RuntimeError`.
- Reindex shared interval endpoints consistently and restore generator index
  maps when loading HDF5 files.
- Handle empty diagrams, infinite intervals, and multi-digit homological degrees
  in plots. Apply persistence thresholds before truncating displayed
  coordinates, honor scatter colour overrides, and fix colour argument typing
  and a Matplotlib warning.
- Fix C++ build compatibility with Eigen 5 and GCC in C++23 mode.

### Documentation and development

- Update the example notebook for the persistence API, representatives, natural
  maps, plotting, and HDF5. Upgrade it to notebook format 4.5 with cell IDs.
- Update Sphinx API references, omit notebook checkpoints from API discovery,
  and refresh installation and build instructions.
- Separate documentation dependencies from development dependencies and add
  Makefile targets for installation, tests, dependency updates, lockfiles, and
  distributions. Refresh dependency locks, pybind11 3, ConstrainedMiniball,
  vcpkg baselines, CI actions, and package license metadata.
- Expand tests for persistence modules, degree limits, representatives, natural
  maps, HDF5 validation and round trips, and plotting.
