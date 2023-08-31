.. Chalc documentation master file, created by
   sphinx-quickstart on Thu Aug 24 18:00:30 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
.. role:: underline
   :class: underline

Chalc Documentation
===================
Chalc is a Python package for computing the persistent homology of :underline:`ch`\ romatic :underline:`al`\ pha :underline:`c`\ omplexes of coloured point clouds in Euclidean space. Chalc is written in C++ and relies on the `Computational Geometry Algorithms Library (CGAL) <https://www.cgal.org>`_ for fast and accurate geometric computations, and `Phimaker <https://github.com/tomchaplin/phimaker>`_ for persistent homology calculations. 

Documentation index
-------------------
.. toctree::
   :name: toc
   :maxdepth: 2

   theory
   installation
   example
   API

.. sidebar-links::
   :github: 
   :pypi: chalc

   license


See Also
--------
The Python package `chromatic-tda <https://pypi.org/project/chromatic-tda/>`_ (written by the original authors of :ref:`[1] <di Montesano et. al.>`) provides similar functionality, but is currently limited to points in two dimensions with upto three colours. 