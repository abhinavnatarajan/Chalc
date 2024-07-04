# Chalc
![PyPI - Version](https://img.shields.io/pypi/v/chalc?labelColor=222222)
![Documentation](https://img.shields.io/badge/docs-stable-blue?labelColor=222222&link=https%3A%2F%2Fabhinavnatarajan.github.io%2FChalc)
![GitHub Workflow Status (with event)](https://img.shields.io/github/actions/workflow/status/abhinavnatarajan/Chalc/build.yml?labelColor=222222)

Chalc is a Python package to compute chromatic Delaunay filtrations from labelled datasets in Euclidean space.
Chromatic Delaunay filtrations are a family of combinatorial objects that capture spatial relations among the classes of a labelled dataset.
You can find more details on the math behind chromatic Delaunay filtrations in the following papers:

[Chromatic Alpha Complexes](https://arxiv.org/abs/2212.03128)\
[Morse Theory for Chromatic Delaunay Triangulations](https://arxiv.org/abs/2405.19303)


For instructions on how to install and use chalc, read the [documentation](https://abhinavnatarajan.github.io/Chalc).

> [!IMPORTANT]
> Since the release of NumPy 2.0, versions of chalc older than 2.0.0 won't work correctly. This is due to a bug in these package versions related to dependency specification. If you are using an older version of chalc, please make sure to install a compatible version of NumPy (≥ 1.24.2 and < 2.0.0).
