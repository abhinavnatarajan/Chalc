Installation
============

Using pip
---------
The recommended way to download and install chalc is from the `PyPI repository <https://pypi.org/project/chalc/>`_ using pip. Pre-packaged binary distributions are available for Windows and Linux (x86_64), and macOS (x86_64 and aarch64). Chalc works with CPython and PyPy ≥ 3.12.

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

The recommended way to obtain and manage these dependencies is using vcpkg (see the build dependencies section) . It is also recommended to use a `Python virtual environment <https://docs.python.org/3/tutorial/venv.html>`_ for the build process to avoid polluting the package namespace.

Build dependencies
^^^^^^^^^^^^^^^^^^
1. CMake (version 3.23 or later).
2. On Windows: Visual Studio 2019 or later.
    On Linux: GCC11 or later.
    On MacOS: Clang 14 or later.
3. (Recommended) `Microsoft vcpkg <https://vcpkg.io/>`_ C++ dependency manager.
4. (MacOS only) The build tools automake, autoconf, and libtool. You can install these with ``brew install automake autoconf libtool``.

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

3. Build the package using your build tool of choice.

.. tab-set::

    .. tab-item:: uv
        :sync: uv

        .. code-block:: bash

            # To build chalc as an editable package in a new virtual environment.
            # This will also install the dependencies for testing and building documentation.
            uv sync --locked --all-groups --exact --verbose

    .. tab-item:: pip
        :sync: pip

        .. code-block:: bash

            # To install chalc into the current Python environment from the package lock file.
            # This will also install the dependencies for testing and building documentation.
            pip lock -r pylock.toml

4. (Optional) If you have GNU Make installed, you can generate type stubs for the package by running ``make stubs`` in the package root directory. If you do not have GNU Make installed, you can run the following commands:

   .. code-block:: bash

            python -m pybind11_stubgen chalc.chromatic --numpy-array-use-type-var --output-dir ./src
            python -m pybind11_stubgen chalc.filtration --numpy-array-use-type-var --output-dir ./src

Building the Documentation
--------------------------

To build the documentation, the development dependencies of the project need to be installed into the current environment.
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
            python -m pybind11_stubgen chalc.chromatic --numpy-array-use-type-var --output-dir ..\src
            python -m pybind11_stubgen chalc.filtration --numpy-array-use-type-var --output-dir ..\src
            sphinx-build -M html source build

This will build the documentation into the folder ``docs/build`` with root ``index.html``.

