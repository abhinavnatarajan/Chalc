# Changelog

## 15.0.1

### Fixed
- Prevent the legend in `plot_sixpack` from overlapping with the axes of diagrams.

## 15.0.0

Changes since 14.2.0.

### Breaking API changes

- Move persistence types and filtration morphisms from `chalc.sixpack` to
  `chalc.morphisms`.
- Rename `FiltrationQuotient` to `FiltrationGluing`, `SubChromaticQuotient` to
  `SubChromaticGluing`, and `KChromaticQuotient` to `KChromaticGluing`.
- Replace `SimplexPairings` with `AugmentedPersistenceModule`, which associates
  interval pairings with endpoint entrance times and homological degrees. Use
  `module.intervals(degrees)` to obtain birth/death arrays: an integer selects
  one array, while a sequence or `None` returns a dictionary indexed by degree.
  `module.pairings` is a read-only mapping from birth indices to death indices,
  with `None` representing infinite death. `module.degrees` is a read-only
  mapping from birth indices to homological degrees.
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
- Store representative-cycle supports and natural-map columns in
  `SixPackRepresentatives` and `SixPackMorphisms` as one-dimensional `np.int64`
  arrays. Relative representatives are pairs of these arrays. Direct
  construction should supply array values, for example
  `np.asarray(indices, dtype=np.int64)`.
- Return an `np.int64` array from `AugmentedPersistenceModule.nonempty_degrees()`
  instead of a list.
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
  instead of `dimensions`. Replace `draw_filtration`'s `time` and
  `include_colours` arguments with an optional `select(simplex, index)` callback,
  where `index` is the simplex's zero-based boundary-matrix index. The remaining
  plotting options are keyword-only, as are `animation_length` and subsequent
  arguments to `animate_filtration`.
- Require `phimaker >= 0.5.0`.

### Added

- Export `IndexArray` as the shared type for one-dimensional `np.int64` arrays.
  Accept NumPy degree arrays and integer scalars in interval queries and plotting,
  and NumPy integer scalars when looking up an interval by birth index.
- Add `sixpack(get_representatives=True)` to compute representative cycles for
  all six modules and the natural maps between kernel, domain, image, codomain,
  and cokernel. Access them through `pack.representatives` and `pack.morphisms`.
  Generator indices are local to each module. For inclusions, all representative
  supports use codomain boundary-matrix indices. For gluing maps, kernel and
  domain supports use domain indices, while image, codomain, and cokernel
  supports use codomain indices.
- Represent relative cycles by pairs of chains `(a, b)` satisfying `da = 0`
  and `f(a) = db`. For gluing maps, `a` uses domain indices and `b` uses codomain
  indices; for inclusions, both use codomain indices.
- Preserve modules, original pairings, representatives, maps, and the endpoint
  relabelling through HDF5 round trips, including both chains of each relative
  representative. Add `AugmentedPersistenceModule.to_hdf5` for writing individual
  modules. Consumer `pack.metadata` is assigned directly and is not serialized.
  Gluing metadata records the maps into the mapping cylinder, the subfiltration
  containing each domain simplex, and the domain-to-codomain simplex map.
- Add `ChromaticSheaf.from_chromatic_points` to construct persistence modules
  indexed by subsets of colours, restricted by `max_num_colours`. Query modules
  with `stalk` and inclusion maps with `morphism`. Include the zero module at
  the empty colour set, support nonconsecutive colour labels, and cache composed
  maps and identities.
- Add `PersistenceModuleMorphism` for sparse maps over F2 between interval bases,
  with identity and zero maps and composition using `@`. Validate generator
  indices, homological degrees, and interval compatibility; remove coefficients
  between disjoint intervals after composition.
- Add `SubfiltrationInclusion(subfiltration, filtration)` for filtrations with
  the same vertex set. Reject different vertex universes and support
  subfiltrations of lower dimension. Shared simplices must have matching
  filtration times.
- Include all boundary faces of simplices selected by `draw_filtration` by
  default. Pass `include_faces=False` to draw only the selected simplices.
- Expose separate feature and truncation legend controls in
  `plot_persistence_diagram`, which consistently returns its axes, including for
  empty diagrams. `plot_sixpack` calculates default truncations for module names
  omitted from a supplied mapping.
- Declare free-threaded Python support in both native extensions. Independent
  computations and concurrent reads of unchanged data can run without caller
  locks; conflicting reads and writes to shared data require a common lock.

### Fixed

- Apply homological degree limits correctly to kernel intervals and the relative
  homology of gluing maps, accounting for shifted degrees in the mapping cylinder.
- Validate interval indices, entrance-time ordering, and degrees before
  filtering; reject representatives and maps that refer to invalid generators.
  Reject malformed HDF5 datasets with a `RuntimeError`, including mismatched
  generator keys in the two components of relative representatives.
- Reindex shared interval endpoints consistently and restore generator index
  maps when loading HDF5 files.
- Handle empty diagrams, infinite intervals, and multi-digit homological degrees
  in plots. Apply persistence thresholds before truncating displayed
  coordinates, honor scatter colour overrides, and fix colour argument typing
  and a Matplotlib warning.
- Draw each selected simplex from its own boundary-matrix label, preventing
  unrelated triangles from appearing when drawing selected cycles and faces.
- Solve induced persistence maps over F2 with empty or dependent columns and
  unsorted supports; reject inconsistent systems and discard coefficients
  between disjoint intervals after solving.
- Fix C++ build compatibility with Eigen 5 and GCC in C++23 mode.

### Array storage and processing

- Keep endpoint-index helpers, degree selection, and sparse-solver results as
  NumPy arrays. Retain lists for incrementally built Rust inputs and sets for
  mutable F2 elimination work.
- Write sparse HDF5 supports using array operations and load rows as array views
  into the flat support buffers. Store the two components of relative
  representatives in separate sparse groups.
- Compare array-valued persistence containers by their contents, preserving
  equality checks and HDF5 round-trip comparisons.

### Documentation and development

- Make GUDHI comparisons optional through the `test-gudhi` dependency group.
  The standard test group and development installs no longer require GUDHI;
  only the five marked comparison tests skip when it is absent. Add
  `make tests-gudhi` to run the complete suite with GUDHI installed.
- Update the example notebook for the persistence API, representatives, natural
  maps, plotting, and HDF5, including examples of drawing selected representative
  cycles. Upgrade it to notebook format 4.5 with cell IDs.
- Add a thread-safety documentation page covering independent computations,
  shared reads, caller locks, plotting, HDF5, and free-threaded Python. Test
  concurrent computations, reads, and synchronized updates, and check that
  importing chalc leaves the GIL disabled on free-threaded builds.
- Update Sphinx API references, omit notebook checkpoints from API discovery,
  and refresh installation and build instructions.
- Separate documentation dependencies from development dependencies and add
  Makefile targets for installation, tests, dependency updates, lockfiles, and
  distributions. Refresh dependency locks, pybind11 3, ConstrainedMiniball,
  vcpkg baselines, CI actions, and package license metadata.
- Add Lefthook pre-commit formatting for Python, C++, TOML, and notebook output
  cleanup, with explicit dependency groups for the hook tools.
- Expand tests for persistence modules, degree limits, representatives, natural
  maps, HDF5 validation and round trips, and plotting. Cover sheaf construction,
  induced maps, composition, interval validation, and selected-simplex geometry.
