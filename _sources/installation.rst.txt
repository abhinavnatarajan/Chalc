Installation
============


Using pip
---------
The recommended way to download and install chalc is from the `PyPI repository <https://pypi.org/project/chalc/>`_ using pip. Prepackaged binary distributions are available for 64-bit Windows, macOSX, and Linux. Chalc requires CPython or PyPy \>=3.8.

.. code-block:: console

   $ pip install chalc

Build from source
-----------------
If a pre-packaged binary distribution is not available for your platform, you can build chalc from source.

Dependencies
^^^^^^^^^^^^
Chalc is a C++ extension module for Python and has several additional dependencies.

1. The `Eigen C++ library <https://eigen.tuxfamily.org/index.php?title=Main_Page>`_ (tested with version 3.4.0).
2. The `GNU MP Library <https://gmplib.org/>`_ (tested with version 6.3.1) and the `GNU MPFR Library <https://www.mpfr.org/>`_ (tested with version 4.2.0) for exact geometric computation. 
3. The `Computational Geometry Algorithms Library (CGAL) <https://www.cgal.org/>`_ library (version 5.6).
4. The `Boost C++ libraries <https://www.boost.org/>`_ (transitive dependency through CGAL). 
   
The recommended way to obtain and manage these dependencies is using the `vcpkg <https://vcpkg.io/>`_ C++ dependency manager. It is also recommended to use a `Python virtual environment <https://docs.python.org/3/tutorial/venv.html>`_ for the build process to avoid polluting the package namespace. 

Build steps
^^^^^^^^^^^

Clone the git repository. 

.. code-block:: console

   $ git clone https://github.com/abhinavnatarajan/chalc
   $ cd chalc

If you have vcpkg installed on your system, make sure that the environment variable ``VCPKG_INSTALLATION_ROOT`` is set and points to the base directory where vcpkg is installed.
   
.. tab-set::

   .. tab-item:: Bash
      :sync: bash

      .. code-block:: bash

         $ export VCPKG_INSTALLATION_ROOT=/path/to/vcpkg/dir

   .. tab-item:: Powershell
      :sync: powershell

      .. code-block:: powershell

         > $Env:VCPKG_INSTALLATION_ROOT = 'C:\path\to\vcpkg\dir'

If you do not have vcpkg installed, the build process will automatically download vcpkg into a temporary directory and fetch the required dependencies.

.. note:: 
   If you would like to disable the use of vcpkg altogether, set the environment variable ``NO_USE_VPKG``. You will have to ensure that all build requirements are met and the appropriate entries are recorded in ``PATH`` (see the file `CMakeLists.txt <https://github.com/abhinavnatarajan/Chalc/blob/master/CMakeLists.txt>`_ for details).
   
   .. tab-set::

      .. tab-item:: Bash
         :sync: bash

         .. code-block:: bash

            $ export NO_USE_VPKG
      
      .. tab-item:: Powershell
         :sync: powershell

         .. code-block:: powershell

            > $Env:NO_USE_VPKG = $null

Build the package wheel using ``pip wheel`` and install the compiled binary.

.. code-block:: console

   $ pip wheel . -w outputdir
   $ pip install outputdir/<name_of_generated_wheel>.whl