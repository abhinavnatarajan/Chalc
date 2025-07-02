..  Chalc documentation master file, created by
    sphinx-quickstart on Thu Aug 24 18:00:30 2023.
    You can adapt this file completely to your liking, but it should at least
    contain the root `toctree` directive.
.. role:: underline
    :class: underline

Chalc Documentation
===================
Chalc is a Python package for computing chromatic Delaunay filtrations of labelled point clouds in Euclidean space, and associated 6-packs of persistent homology diagrams.
6-packs are algebraic structures which can be used to quantify spatial relationships between various classes of points in a labelled point cloud, and chromatic Delaunay filtrations are a computational tool that make it feasible to compute 6-packs for low-dimensional labelled datasets.
Chalc is named for :underline:`ch`\ romatic :underline:`al`\ pha :underline:`c`\ omplexes, but it can also compute chromatic Delaunay--Čech and chromatic Delaunay--Rips filtrations.
Chalc is written in C++ and relies on the `Computational Geometry Algorithms Library (CGAL) <https://www.cgal.org>`_ for geometric computations, `Phimaker <https://github.com/tomchaplin/phimaker>`_ for persistent homology calculations, and `Intel OneAPI Threading Building Blocks (TBB) <https://www.threadingbuildingblocks.org>`_ for parallelism.

Here are some key features of Chalc:

#. Internal computations are performed with multiple precision exact arithmetic.
#. Chalc will try to produce reasonable results even when the input points are degenerate.
#. Parallelism is used wherever possible to speed up computations (unless explicitly disabled).
#. Lock-free computation of persistence diagrams using the algorithm of :ref:`[2] <2>`, along with the clearing optimisation from :ref:`[3] <3>`.

There are no constraints on the dimensionality of the point clouds that Chalc can handle, and up to 16 colours are supported.
However, the package is primarily intended for 2D and 3D point clouds with two or three colours and a few thousand points.
You may run into significant memory bottlenecks in higher dimensions or with more colours unless you have very few points.

Contributing
------------
Writing a mathematical library that is performant, bug-free, composable, and well-documented is a considerable task, and we welcome contributions from users of this package.
If you find a bug, or an error in the documentation, or have a feature request, or would like to contribute, please open an `issue on the GitHub repository <https://github.com/abhinavnatarajan/Chalc/issues>`_.
If you have a question about how to use the package, please open a `discussion on the GitHub repository <https://github.com/abhinavnatarajan/Chalc/discussions>`_.

See Also
--------
The Python package `chromatic-tda <https://pypi.org/project/chromatic-tda/>`_ (written by the authors of :ref:`[1] <1>`) has similar functionality, and can compute chromatic alpha filtrations and their six-packs of persistence diagrams, but does not provide functionality to compute chromatic Delaunay--Čech or chromatic Delaunay--Rips filtrations.

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


References
----------
#. .. _1:
   Sebastiano Cultrera di Montesano, Ondřej Draganov, Herbert Edelsbrunner, and Morteza Saghafian. "Chromatic Alpha Complexes". Feb. 7, 2024. arXiv: `2212.03128 <https://arxiv.org/abs/2212.03128>`_.
#. .. _2:
   Dmitriy Morozov and Arnur Nigmetov. 2020. Towards Lockfree Persistent Homology. In Proceedings of the 32nd ACM Symposium on Parallelism in Algorithms and Architectures (SPAA '20). Association for Computing Machinery, New York, NY, USA, 555–557. https://doi.org/10.1145/3350755.3400244.
#. .. _3:
   Bauer, U., Kerber, M., Reininghaus, J. (2014). Clear and Compress: Computing Persistent Homology in Chunks. In: Bremer, PT., Hotz, I., Pascucci, V., Peikert, R. (eds) Topological Methods in Data Analysis and Visualization III. Mathematics and Visualization. Springer, Cham. https://doi.org/10.1007/978-3-319-04099-8_7
