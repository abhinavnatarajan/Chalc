# Build dependencies, must be listed in pyproject.toml
from skbuild import setup
from setuptools import find_packages
from platform import python_version
import re

# Available through Python standard library
import json
from pathlib import Path
from shutil import rmtree

PROJECT_SOURCE_DIR = Path(__file__).parent
# For some reason, running this file twice in a row causes the build to fail:
# fatal error C1083: Cannot open include file: 'io.h'
# Therefore the workaround is to clean the `_skbuild` directory before running
SKBUILD_DIR = PROJECT_SOURCE_DIR / "_skbuild"
if SKBUILD_DIR.exists():
    print(f"Removing previous installation: {SKBUILD_DIR}")
    rmtree(str(SKBUILD_DIR))

# Get vcpkg commit id to pass to CMake
with open(PROJECT_SOURCE_DIR / "vcpkg.json") as f:
    vcpkg_json = json.load(f)
    # Required
    VCPKG_COMMIT_ID = vcpkg_json["builtin-baseline"]

# Get name and version for CMake from pyproject.toml
if python_version() >= '3.11.0':
    from tomllib import load as tomlread
    with open(PROJECT_SOURCE_DIR / "pyproject.toml", "rb") as f:
        proj_props = tomlread(f)['project']
else:
    from toml import loads as tomlread
    with open(PROJECT_SOURCE_DIR / "pyproject.toml") as f:
        proj_props = tomlread(f.read())['project']

project = proj_props['name']
# Convert PEP440 version info to SemVer
def pep440toSemVer(x):
    res = []
    for v in x.split('.'):
        if (re.search('post|pre|dev', v) is not None):
            break
        res.append(re.sub('a.+|b.+|rc.+', '', v))
    return '.'.join(res[0:3])

release = pep440toSemVer(proj_props['version'])

# scikit-build will take care of puting our compiled C++ library together with
# our python package so it can access it. The name of the python package will
# be determined by the name of the folder that contains an `__init__.py` file.
# In this repository, python packages must be placed under path defined by
# `python_packages_root`.
# In order to change the name of the package, the name of the folder that
# contains the `__init__.py` file must be changed.
python_packages_root = "src/python"
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