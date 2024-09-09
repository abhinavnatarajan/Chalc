.. Chalc documentation master file, created by
   sphinx-quickstart on Thu Aug 24 18:00:30 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
.. role:: underline
   :class: underline

Chalc Documentation
===================
Chalc is a Python package for computing chromatic Delaunay filtrations of labelled point clouds in Euclidean space, and their persistent homology.
Chalc is named for :underline:`ch`\ romatic :underline:`al`\ pha :underline:`c`\ omplexes, but it can also compute chromatic Delaunay--Čech and chromatic Delaunay--Rips filtrations.
Chalc is written in C++ and relies on the `Computational Geometry Algorithms Library (CGAL) <https://www.cgal.org>`_ for fast and accurate geometric computations, and `Phimaker <https://github.com/tomchaplin/phimaker>`_ for persistent homology calculations.

Documentation index
-------------------
.. toctree::
   :name: toc
   :maxdepth: 2

   installation
   example
   API
   license

.. sidebar-links::
   :github:
   :pypi: chalc

Citing this Software
--------------------
If you want to use chalc in your work, please use the following citation:

.. tab-set::

   .. tab-item:: Text
      :sync: text

      Abhinav Natarajan, Thomas Chaplin, Adam Brown, and Maria-Jose Jimenez. "Morse Theory for Chromatic Delaunay Triangulations". May 30, 2024. arXiv: `2405.19303 <https://arxiv.org/abs/2405.19303>`_.

   .. tab-item:: BibTex
      :sync: bibtex

      .. code-block:: bibtex

         @misc{natarajan2024morse,
            title={Morse Theory for Chromatic Delaunay Triangulations},
            author={Abhinav Natarajan and Thomas Chaplin and Adam Brown and Maria-Jose Jimenez},
            year={2024},
            eprint={2405.19303},
            archivePrefix={arXiv},
            primaryClass={math.AT}
         }

See Also
--------
The Python package `chromatic-tda <https://pypi.org/project/chromatic-tda/>`_ (written by the original authors of :ref:`[1] <Sebastiano>`) provides similar functionality, but is currently limited to points in two dimensions with upto three colours, although an alpha version with support for higher dimensions is available.

References
----------

#. _`Sebastiano` Cultrera di Montesano, Ondřej Draganov, Herbert Edelsbrunner, and Morteza Saghafian. "Chromatic Alpha Complexes". Feb. 7, 2024. arXiv: `2212.03128 <https://arxiv.org/abs/2212.03128>`_.
