.. Chalc documentation master file, created by
   sphinx-quickstart on Thu Aug 24 18:00:30 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
.. role:: underline
   :class: underline

Chalc Documentation
===================
Chalc is a Python package for computing the persistent homology of :underline:`ch`\ romatic :underline:`al`\ pha :underline:`c`\ omplexes of point clouds in Euclidean space. Chalc is written in C++ and relies on the `Computational Geometry Algorithms Library (CGAL) <https://www.cgal.org>`_ for fast and accurate geometric computations, and the `Lockfree Persistent Homology Algorithm Toolbox (loPHAT) <https://github.com/tomchaplin/lophat>`_ for persistent homology calculations. 

Documentation index
-------------------
.. toctree::
   :name: toc
   :maxdepth: 2

   installation
   example
   API

.. sidebar-links::
   :github: 
   :pypi: chalc

   license

See Also
--------
The Python package `chromatic-tda <https://pypi.org/project/chromatic-tda/>`_ provides similar functionality, but is currently limited to points in two dimensions with upto three colours.

References
----------

#. di Montesano et. al., “Persistent Homology of Chromatic Alpha Complexes”.
   Online preprint available at `<https://arxiv.org/abs/2212.03128>`_.
   Accessed: 2023-02-28 22:07:23 UTC.
   DOI: 10.48550/arXiv.2212.03128.

#. E. Welzl, “Smallest enclosing disks (balls and ellipsoids),” 
   in New Results and New Trends in Computer Science, H. Maurer, Ed., 
   in Lecture Notes in Computer Science. Berlin, Heidelberg: Springer, 
   1991, pp. 359–370. doi: `10.1007/BFb0038202 <https://doi.org/10.1007/BFb0038202>`_.