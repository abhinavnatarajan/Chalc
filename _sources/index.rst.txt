..  Chalc documentation master file, created by
    sphinx-quickstart on Thu Aug 24 18:00:30 2023.
    You can adapt this file completely to your liking, but it should at least
    contain the root `toctree` directive.
.. role:: underline
    :class: underline

Chalc Documentation
===================
Chalc is a Python package for computing chromatic Delaunay filtrations of labelled point clouds in Euclidean space, and the associated 6-packs of persistent homology diagrams.
6-packs are algebraic structures which can be used to quantify spatial relationships between various classes of points in a labelled point cloud, and chromatic Delaunay filtrations are a computational tool that make it feasible to compute 6-packs for low-dimensional labelled datasets.
Chalc is named for :underline:`ch`\ romatic :underline:`al`\ pha :underline:`c`\ omplexes, which were defined in :ref:`[1] <1>`, but the package can also compute chromatic Delaunay--Čech and chromatic Delaunay--Rips filtrations, which are computationally simpler and are based on the ideas in :ref:`[2] <2>`.
Chalc is written in C++ and relies on the `Computational Geometry Algorithms Library (CGAL) <https://www.cgal.org>`_ for geometric computations, `Phimaker <https://github.com/tomchaplin/phimaker>`_ for persistent homology calculations, and `Intel OneAPI Threading Building Blocks (TBB) <https://www.threadingbuildingblocks.org>`_ for parallelism.

Here are some key features of Chalc:

#. Internal computations are performed with multiple precision exact arithmetic.
#. Chalc will try to produce reasonable results even when the input points are degenerate.
#. Parallelism is used wherever possible to speed up computations (unless explicitly disabled).
#. Lock-free computation of persistence diagrams using the algorithm of :ref:`[3] <3>`, along with the clearing optimisation from :ref:`[4] <4>`.

There are no constraints on the dimensionality of the point clouds that Chalc can handle, and up to 16 colours are supported.
However, the package is primarily intended for 2D and 3D point clouds with two or three colours and a few thousand points.
You may run into significant memory bottlenecks in higher dimensions or with more colours unless you have very few points.

Contributing
------------
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

.. _Morse Theory for Chromatic Delaunay Triangulations:
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
The Python package `chromatic-tda <https://pypi.org/project/chromatic-tda/>`_ is written and maintained by the authors of :ref:`[1] <1>`, the original creators of chromatic alpha filtrations, and we are indebted to them for their work.
Chromatic-tda can compute chromatic alpha filtrations and their six-packs of persistence diagrams.
It can also compute cycle representatives for persistent features, which chalc does not currently provide (although this is planned for a future release).

On the other hand, chalc can compute chromatic Delaunay--Čech or chromatic Delaunay--Rips filtrations, which chromatic-tda does not currently support.
Chalc is generally faster than chromatic-tda due to its lower-level implementation and the use of parallelism, and you may find this useful if you are working with a large number of datasets or datasets with a large number of points.
The following is a chart comparing the runtime performance of chalc v\ |chalc-version| and chromatic-tda v\ |chromatic-tda-version| on randomly generated 2-dimensional point clouds with two and three colours.
These tests were run on a laptop with an AMD Ryzen 7 6800U processor (8 cores @2.7Ghz, x86-64 ISA) and 16 GB of RAM, running Ubuntu 24.04 (kernel 6.12.10-76061203-generic) and using Python 3.13.2.
The code used to run these tests and generate the charts are available in the ``docs/source`` folder of the chalc repository, and the raw results are available as a JSON file in the same folder.

.. image:: /timing_results.png


References
----------
#. .. _1:

   Sebastiano Cultrera di Montesano, Ondřej Draganov, Herbert Edelsbrunner, and Morteza Saghafian. "Chromatic Alpha Complexes". Feb. 7, 2024. arXiv: `2212.03128 <https://arxiv.org/abs/2212.03128>`_.

#. .. _2:

   Abhinav Natarajan, Thomas Chaplin, Adam Brown, and Maria-Jose Jimenez. "Morse Theory for Chromatic Delaunay Triangulations". May 30, 2024. arXiv: `2405.19303 <https://arxiv.org/abs/2405.19303>`_.

#. .. _3:

   Dmitriy Morozov and Arnur Nigmetov. 2020. Towards Lockfree Persistent Homology. In Proceedings of the 32nd ACM Symposium on Parallelism in Algorithms and Architectures (SPAA '20). Association for Computing Machinery, New York, NY, USA, 555–557. `<https://doi.org/10.1145/3350755.3400244>`_.

#. .. _4:

   Bauer, U., Kerber, M., Reininghaus, J. (2014). Clear and Compress: Computing Persistent Homology in Chunks. In: Bremer, PT., Hotz, I., Pascucci, V., Peikert, R. (eds) Topological Methods in Data Analysis and Visualization III. Mathematics and Visualization. Springer, Cham. `<https://doi.org/10.1007/978-3-319-04099-8_7>`_.
