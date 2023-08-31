# Build dependencies, must be listed in pyproject.toml
from skbuild import setup
from setuptools import find_packages
from platform import python_version
from packaging.version import parse

# Available through Python standard library
import json
from pathlib import Path
from shutil import rmtree

project_root_dir = Path(__file__).parent
# For some reason, running this file twice in a row causes the build to fail:
# fatal error C1083: Cannot open include file: 'io.h'
# Therefore the workaround is to clean the `_skbuild` directory before running
SKBUILD_DIR = project_root_dir / "_skbuild"
if SKBUILD_DIR.exists():
    print(f"Removing previous installation: {SKBUILD_DIR}")
    rmtree(str(SKBUILD_DIR))

# Get vcpkg commit id to pass to CMake
with open(project_root_dir / "vcpkg.json") as f:
    vcpkg_json = json.load(f)
    # Required
    VCPKG_COMMIT_ID = vcpkg_json["builtin-baseline"]

# Get name and version for CMake from pyproject.toml
if parse(python_version()) >= parse('3.11'):
    from tomllib import loads as tomlread
else:
    from toml import loads as tomlread
with open(project_root_dir / "pyproject.toml") as f:
    proj_props = tomlread(f.read())['project']

project = proj_props['name']
release = '.'.join(map(str, parse(proj_props['version']).release))

# scikit-build will take care of puting our compiled C++ library together with
# our python package so it can access it. The name of the python package will
# be determined by the name of the folder that contains an `__init__.py` file.
# In this repository, python packages must be placed under path defined by
# `python_packages_root`.
# In order to change the name of the package, the name of the folder that
# contains the `__init__.py` file must be changed.
python_packages_root = project_root_dir / 'src' / 'python'
packages = find_packages(python_packages_root)

setup(
    packages = packages,
    package_dir = {"" : python_packages_root},
    cmake_install_dir = python_packages_root + "/" + packages[0],
    cmake_with_sdist = False,
    package_data = { packages[0] : ["*.dll"] },
    cmake_args=["-DSKBUILD_PROJECT_NAME:STRING=" + project,
                "-DSKBUILD_PROJECT_VERSION:STRING=" + release,
                "-DVCPKG_COMMIT_ID:STRING=" + VCPKG_COMMIT_ID,
                "-DCMAKE_BUILD_TYPE:STRING=Release"]
)