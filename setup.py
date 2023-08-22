# Build dependencies, must be listed in pyproject.toml
from skbuild import setup
from setuptools import find_packages

# Available through Python standard library
import json
from pathlib import Path
from shutil import rmtree
from os import environ
from platform import system

PROJECT_SOURCE_DIR = Path(__file__).parent
# For some reason, running this file twice in a row causes the build to fail:
# fatal error C1083: Cannot open include file: 'io.h'
# Therefore the workaround is to clean the `_skbuild` directory before running
SKBUILD_DIR = PROJECT_SOURCE_DIR / "_skbuild"
if SKBUILD_DIR.exists():
    print(f"Removing previous installation: {SKBUILD_DIR}")
    rmtree(str(SKBUILD_DIR))

# In order to avoid specifying package name and version in multiple files, we
# will use `vcpkg.json` in the repository root as reference and extract the
# apropiate variables from there.
with open(PROJECT_SOURCE_DIR / "vcpkg.json") as f:
    vcpkg_json = json.load(f)
    # Required
    PROJECT_VERSION_STRING = vcpkg_json["version-semver"]
    PROJECT_NAME = vcpkg_json["name"]
    VCPKG_COMMIT_ID = vcpkg_json["builtin-baseline"]

# For a CI non-Windows build we will change the GMP requirement to "fat"
# This passes the --enable-fat flag while configuring GMP
# See https://gmplib.org/manual/Build-Options
# This ensures that wheels are built for a variety of platforms
if environ.get("CI") and system() != 'Windows':
    for dep in vcpkg_json["dependencies"]:
        if dep["name"] == "gmp":
            dep["features"] = ["fat"]
            break
    with open(PROJECT_SOURCE_DIR / "vcpkg.json", "w") as f:
        json.dump(vcpkg_json, f, indent = 4)

# Read the project description from README.md
long_description = (PROJECT_SOURCE_DIR / "README.md").read_text()
    
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
    name = PROJECT_NAME,
    version = PROJECT_VERSION_STRING,
    long_description = long_description,
    long_description_content_type = "text/markdown",
    author = "Abhinav Natarajan",
    url = "https://github.com/abhinavnatarajan/chalc",
    download_url = "https://github.com/abhinavnatarajan/chalc",
    license = "GPL",
    packages = packages,
    package_dir = {"" : python_packages_root},
    cmake_install_dir = python_packages_root + "/" + packages[0],
    cmake_with_sdist = False,
    package_data = { packages[0] : ["*.dll"] },
    python_requires = ">=3.8",
    cmake_args=["-DSKBUILD_PROJECT_NAME:STRING=" + PROJECT_NAME,
                "-DSKBUILD_PROJECT_VERSION:STRING=" + PROJECT_VERSION_STRING,
                "-DVCPKG_COMMIT_ID:STRING=" + VCPKG_COMMIT_ID,
                "-DCMAKE_BUILD_TYPE:STRING=Release"]
)