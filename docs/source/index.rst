..  Chalc documentation master file, created by
    sphinx-quickstart on Thu Aug 24 18:00:30 2023.
    You can adapt this file completely to your liking, but it should at least
    contain the root `toctree` directive.
.. role:: underline
    :class: underline

Chalc Documentation
===================
Chalc is a Python package for computing chromatic Delaunay filtrations of labelled point clouds in Euclidean space, and associated six-packs of persistent homology diagrams.
Chalc is named for :underline:`ch`\ romatic :underline:`al`\ pha :underline:`c`\ omplexes, but it can also compute chromatic Delaunay--Čech and chromatic Delaunay--Rips filtrations.
Chalc is written in C++ and relies on the `Computational Geometry Algorithms Library (CGAL) <https://www.cgal.org>`_ for geometric computations, `Phimaker <https://github.com/tomchaplin/phimaker>`_ for persistent homology calculations, and `Intel OneAPI Threading Building Blocks (TBB) <https://www.threadingbuildingblocks.org>`_ for parallelism.

Chalc uses exact arithmetic internally to compute filtration values and will produce reasonable results even for point clouds with degeneracies.
Chalc also uses parallelism wherever possible to speed up computations unless parallelism is explicitly disabled.

Contributing
------------
Chalc is a work in progress.
Writing a mathematical library that is performant, bug-free, composable, and well-documented is a considerable task, and we welcome contributions from users of this package.
If you find a bug, or an error in the documentation, or have a feature request, or would like to contribute, please open an `issue on the GitHub repository <https://github.com/abhinavnatarajan/Chalc/issues>`_.
If you have a question about how to use the package, please open a `discussion on the GitHub repository <https://github.com/abhinavnatarajan/Chalc/discussions>`_.

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
If you find chalc useful, please consider citing it as below:

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
The Python package `chromatic-tda <https://pypi.org/project/chromatic-tda/>`_ (written by the original authors of :ref:`[1] <Sebastiano>`) has similar functionality, and can compute chromatic alpha filtrations and their six-packs of persistence diagrams, but does not provide functionality to compute chromatic Delaunay--Čech or chromatic Delaunay--Rips filtrations.

References
----------

#. _`Sebastiano` Cultrera di Montesano, Ondřej Draganov, Herbert Edelsbrunner, and Morteza Saghafian. "Chromatic Alpha Complexes". Feb. 7, 2024. arXiv: `2212.03128 <https://arxiv.org/abs/2212.03128>`_.

