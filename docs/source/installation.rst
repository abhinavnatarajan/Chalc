Installation
============

Using pip
---------
The recommended way to download and install chalc is from the `PyPI repository <https://pypi.org/project/chalc/>`_ using pip. Pre-packaged binary distributions are available for Windows and Linux with x86_64 CPUs, and MacOS with M-series CPUs.

.. tab-set::

    .. tab-item:: uv
        :sync: uv

        .. code-block:: bash

            # To chalc to your uv project
            uv add chalc

    .. tab-item:: pip
        :sync: pip

        .. code-block:: bash

            # To install chalc into the current Python environment
            pip install chalc


Build from source
-----------------
If a pre-packaged binary distribution is not available for your platform, you can build chalc from source.

Dependencies
^^^^^^^^^^^^
Chalc is a C++ extension module for Python and has several additional dependencies.

1. `Eigen C++ library <https://eigen.tuxfamily.org/index.php?title=Main_Page>`_ (tested with version 3.4.0).
2. `GNU MP Library <https://gmplib.org/>`_ (tested with version 6.3.1) and the `GNU MPFR Library <https://www.mpfr.org/>`_ (tested with version 4.2.0) for exact geometric computation.
3. `Computational Geometry Algorithms Library (CGAL) <https://www.cgal.org/>`_ library (tested with version 6.0.1).
4. `Boost C++ libraries <https://www.boost.org/>`_ (transitive dependency through CGAL).
5. `Intel OneAPI Threading Building Blocks (TBB) <https://www.threadingbuildingblocks.org/>`_ (tested with version 2022.1.0).

The recommended way to obtain and manage these dependencies is using vcpkg (see the build dependencies section).

Build dependencies
^^^^^^^^^^^^^^^^^^
1. CMake (version 4.0.0 or later).
2. On Windows: Visual Studio 2019 or later.
    On Linux: GCC11 or later.
    On MacOS: Clang 14 or later.
3. On Linux and MacOS, the build tools automake, autoconf, autoconf-archive, and libtool. Install them with ``apt install autoconf autoconf-archive automake libtool``, ``dnf install autoconf autoconf-archive automake libtool``, ``pacman -S autoconf autoconf-archive automake libtool``, or ``brew install autoconf autoconf-archive automake libtool``.
4. The Python project tool `uv <https://astral.sh/uv>`_.
5. (Recommended) `Microsoft vcpkg <https://vcpkg.io/>`_ C++ dependency manager.
6. (Recommended) `GNU Make <https://www.gnu.org/software/make/>`_.

Build steps
^^^^^^^^^^^

1. Clone the git repository.

.. code-block:: bash

    git clone https://github.com/abhinavnatarajan/chalc
    cd chalc

2. If you have vcpkg installed on your system, make sure that the environment variable ``VCPKG_ROOT`` is set and points to the base directory where vcpkg is installed.

.. tab-set::

    .. tab-item:: Bash
        :sync: bash

        .. code-block:: bash

            export VCPKG_ROOT=/path/to/vcpkg/dir

    .. tab-item:: Powershell
        :sync: powershell

        .. code-block:: powershell

            $Env:VCPKG_ROOT = 'C:\path\to\vcpkg\dir'

If you do not have vcpkg installed, the build process will automatically download vcpkg into a temporary directory, and fetch and build the required dependencies.

.. note::
    If you would like to disable the use of vcpkg altogether, set the environment variable ``NO_USE_VPKG``. You will have to ensure that all build requirements are met and the appropriate entries are recorded in ``PATH`` (see the file `CMakeLists.txt <https://github.com/abhinavnatarajan/Chalc/blob/master/CMakeLists.txt>`_ for details).

    .. tab-set::

        .. tab-item:: Bash
            :sync: bash

            .. code-block:: bash

                export NO_USE_VPKG

        .. tab-item:: Powershell
            :sync: powershell

            .. code-block:: powershell

                $Env:NO_USE_VPKG = $null

3. Install the project as an editable package, along with the development dependencies.

.. tab-set::

    .. tab-item:: Bash
        :sync: bash

        .. code-block:: bash

            make install

    .. tab-item:: Windows Powershell
        :sync: powershell

        .. code-block:: powershell

            uv sync --verbose --all-groups --no-progress
            uv run python -m pybind11_stubgen chalc.chromatic --numpy-array-use-type-var --output-dir ..\src
            uv run python -m pybind11_stubgen chalc.filtration --numpy-array-use-type-var --output-dir ..\src
            uv lock
            uv export --format pylock.toml --all-groups -o pylock.toml --quiet

Building the Documentation
--------------------------

To build the documentation, the development dependencies of the project need to be installed into the current environment.
You also need to have `GraphViz <https://graphviz.org/download/>`_ installed.
Then run the following commands from the project root directory to build the documentation files.

.. tab-set::

    .. tab-item:: Bash
        :sync: bash

        .. code-block:: bash

            make docs

    .. tab-item:: Windows Powershell
        :sync: powershell

        .. code-block:: powershell

            Set-Location docs
            uv run sphinx-build -M html source build

This will build the documentation into the folder ``docs/build`` with root ``index.html``.
