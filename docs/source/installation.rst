Installation
============

Via pip
-------
The recommended way to install chalc is via pip.
.. code-block:: console

   $ pip install chalc

Build from source
-----------------
Dependencies
^^^^^^^^^^^^
Chalc is a C++ extension module for Python and has several additional dependencies.

1. The `Eigen C++ library <https://eigen.tuxfamily.org/index.php?title=Main_Page>`_ (tested with version 3.4.0).
2. The `GNU MP Library <https://gmplib.org/>`_ (tested with version 6.3.1) and the `GNU MPFR Library <https://www.mpfr.org/>`_ (tested with version 4.2.0) for exact geometric computation. 
3. The `CGAL <https://www.cgal.org/>`_ library (version 5.6).
4. In addition, CGAL depends on the `Boost C++ libraries <https://www.boost.org/>`_. 
   
The recommended way to obtain and manage these dependencies is using the `vcpkg <https://vcpkg.io/>`_ C++ dependency manager. 

Build steps
^^^^^^^^^^^

1. Clone the git repository. 

.. code-block:: console

   $ git clone https://github.com/abhinavnatarajan/chalc
   $ cd chalc

2. If you have vcpkg installed on your system, make sure that the environment variable `VCPKG_INSTALLATION_ROOT` is set and points to the base directory where it is installed.

.. code-block:: console

   $ set




