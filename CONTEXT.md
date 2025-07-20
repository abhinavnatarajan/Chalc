# Chalc Project Workspace

> **Chalc**: Chromatic Delaunay Filtrations for Topological Data Analysis  
> A hybrid Python/C++ library for computing chromatic Delaunay filtrations from labelled datasets in Euclidean space.

## Overview

Chalc is a computational geometry library that combines high-performance C++ algorithms implemented with CGAL, Eigen, Intel OneAPI TBB, and a user-friendly Python interface.
The library focuses on persistent homology and topological data analysis (TDA) through chromatic Delaunay filtrations, enabling spatial relationship analysis in low-dimensional labelled datasets.

**Key Technologies:**
- **C++** for performance-critical computational geometry with CGAL, Eigen, and Intel OneAPI TBB
- **Python** for user-friendly API and visualization
- **pybind11** for seamless Python-C++ integration
- **ConstrainedMiniball** for specialized geometric computations
- **CMake** and **vcpkg** for robust build system
- **GitHub Actions** for automated CI/CD

---

## Functional Groups

### 1. C++ Core Algorithms
*Performance-critical computational geometry engines*

**Focus**: C++ implementation of core algorithms using CGAL, Eigen, and Intel OneAPI TBB for high-performance geometric computations. These components handle the mathematical heavy lifting for chromatic Delaunay triangulations.

**Components:**
- **`src/chalc/chromatic/chromatic.cxx`** - C++ implementation of chromatic alpha complex algorithms using CGAL
- **`src/chalc/chromatic/chromatic.h`** - C++ header defining chromatic alpha complex algorithm interfaces
- **`src/chalc/filtration/filtration.cxx`** - C++ implementation of filtration data type
- **`src/chalc/filtration/filtration.h`** - C++ header defining interfaces for the filtration data type
- **`src/chalc/chromatic_py.cxx`** - Pybind11 bindings exposing C++ chromatic algorithms to Python
- **`src/chalc/filtration_py.cxx`** - Pybind11 bindings exposing C++ filtration type to Python

---

### 2. Python API & Interface
*User-friendly Python layer and integration*

**Focus**: Python interface layer providing accessible access to C++ computational cores. Includes main Python modules, utility functions, and the bridge between Python and C++ through pybind11 bindings.

**Components:**
- **`src/chalc/__init__.py`** - Main Python module initialization and public API exposure
- **`src/chalc/_utils.py`** - Python utility functions and helper methods
- **`src/chalc/chromatic_py.cxx`** - Pybind11 bindings for chromatic algorithms
- **`src/chalc/filtration_py.cxx`** - Pybind11 bindings for filtration algorithms

---

### 3. Geometric Utilities & ConstrainedMiniball
*Specialized geometric computation components*

**Focus**: Geometric utility components, particularly the ConstrainedMiniball library for computing smallest enclosing balls with affine constraints. This is a key computational geometry component supporting the main chromatic algorithms.

**Components:**
- **`src/ConstrainedMiniball/cmb.hpp`** - ConstrainedMiniball library for smallest enclosing balls with affine constraints
- **`src/ConstrainedMiniball/tests/`** - Test suite for ConstrainedMiniball geometric algorithms

---

### 4. Persistent Homology Computations
*Morphisms of persistence modules for topological data analysis*

**Focus**: User-friendly abstractions over morphisms of persistence modules and computation of sixpacks of persistence diagrams. These abstractions are a high-level interface to the low-level API provided by the Python package `phimaker`, which computes the sixpack of persistence diagrams from an annotated boundary matrix. Includes type definitions, morphisms, and mathematical structures for persistent homology computations.

**Components:**
- **`src/chalc/sixpack/__init__.py`** - Sixpack module for mathematical abstractions and type definitions
- **`src/chalc/sixpack/types.py`** - Type definitions and mathematical structures for persistent homology
- **`src/chalc/sixpack/morphisms.py`** - Mathematical morphisms and transformations for topological computations

---

### 5. Visualization & Plotting
*Data visualization and result presentation*

**Focus**: Visualization and plotting capabilities for displaying chromatic Delaunay filtrations, sixpacks of persistent homology diagrams, and geometric structures. Includes matplotlib-based plotting and example visualizations.

**Components:**
- **`src/chalc/plotting.py`** - Matplotlib-based plotting and visualization functions for chromatic complexes

---

### 6. Testing & Validation
*Comprehensive test coverage for quality assurance*

**Focus**: Comprehensive test suite validating both Python and C++ components. Includes unit tests for chromatic computations, filtration algorithms, geometric utilities, and integration tests.

**Components:**
- **`tests/test_chromatic.py`** - Unit tests for chromatic alpha complex computations
- **`tests/test_filtration.py`** - Unit tests for filtration algorithms and computations
- **`tests/test_sixpack.py`** - Unit tests for sixpack mathematical abstractions and types
- **`tests/utils.py`** - Testing utilities and helper functions for test suite
- **`src/ConstrainedMiniball/tests/`** - Test suite for ConstrainedMiniball geometric algorithms

---

### 7. Documentation & Examples
*User guides, API references, and practical examples*

**Focus**: User documentation, API references, mathematical background, and practical examples. Includes Sphinx documentation, Jupyter notebooks demonstrating usage, and comprehensive API documentation.

**Components:**
- **`docs/source/`** - Sphinx documentation source files including API documentation and user guides
- **`docs/source/example.ipynb`** - Comprehensive Jupyter notebook with usage examples
- **`docs/source/API.rst`** - Comprehensive API documentation for all public interfaces
- **`README.md`** - Main project README with overview, installation instructions, and usage examples

---

### 8. Build System & Configuration
*Project setup, building, and deployment infrastructure*

**Focus**: Build system, packaging configuration, and project setup. Includes CMake configuration for C++ components, Python packaging, dependency management, development tools, and CI/CD pipelines.

**Components:**
- **`CMakeLists.txt`** - CMake configuration for building C++ components with CGAL, Eigen, and Intel OneAPI TBB dependencies
- **`pyproject.toml`** - Python packaging configuration with dependencies, build system, and metadata
- **`.`** (root directory) - Development configuration files including .clang-format, .gitignore, and editor configs
- **`vcpkg.json`** - vcpkg dependency configuration for C++ libraries including CGAL and mathematical libraries
- **`.github/`** - GitHub Actions CI/CD workflows for automated testing, building, and deployment

---

## Quick Reference

### Key Directories
```
src/chalc/                    # Main Python package
├── chromatic/               # C++ chromatic algorithms
├── filtration/              # C++ filtration algorithms
├── sixpack/                 # Sixpack computations
├── __init__.py              # Python API
├── _utils.py                # Python utilities
├── plotting.py              # Visualization tools
├── chromatic_py.cxx         # Python bindings (chromatic)
└── filtration_py.cxx        # Python bindings (filtration)

src/ConstrainedMiniball/      # Geometric utilities
tests/                        # Test suite
docs/                         # Documentation
.github/                      # CI/CD workflows
```

### Build System
- **C++**: CMake + vcpkg for dependencies (CGAL, GMP, MPFR, Eigen, TBB)
- **Python**: pyproject.toml with modern packaging
- **CI/CD**: GitHub Actions (build, test, deploy, docs)

### Development Workflow
1. **Core Algorithms**: Implement in C++ with CGAL, Eigen, and TBB
2. **Python Bindings**: Expose via pybind11
3. **Testing**: Comprehensive test coverage
4. **Documentation**: Sphinx + Jupyter notebooks
5. **CI/CD**: Automated testing and deployment
