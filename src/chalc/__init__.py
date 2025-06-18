"""Analyse spatial relationships in low-dimensional data using persistent homology."""

from __future__ import annotations

import platform

from . import filtration

if platform.system() == "Windows":
	# need to load up the GMP and MPFR DLLs
	import os as _os
	from importlib.metadata import files as _files
	from pathlib import Path as _Path

	chalc_files = _files("chalc")
	if chalc_files:
		_bindir = _Path(
			next(
				x for x in chalc_files if len(x.parent.parts) > 0 and x.parent.parts[-1] == "lib"
			).locate(),
		).parent.resolve()
		with _os.add_dll_directory(str(_bindir)):
			from . import chromatic
	else:
		error_msg = (
			"Could not find the required DLLs for chalc on Windows. "
			"Please ensure that the package is installed correctly."
		)
		raise ImportError(error_msg)

else:
	from . import chromatic

from . import plotting, sixpack

__all__ = ["chromatic", "filtration", "plotting", "sixpack"]
