Thread safety
=============

Chalc supports concurrent computations on independent data, and concurrent reads
of shared data that stays unchanged. Chalc does not lock its computational objects
for you. If one thread may change shared data while another reads or changes it,
both threads must use the same lock.
Unsynchronized reads and writes to native filtration data can produce incorrect
results or crash Python.

These rules apply with or without Python's global interpreter lock (GIL). A
free-threaded Python build lets Python threads run at the same time; it does not
make changes to shared objects safe automatically.

Independent computations
------------------------

Each worker can create and change its own inputs and filtration without a lock:

.. code-block:: python

    from concurrent.futures import ThreadPoolExecutor

    import chalc as ch
    import numpy as np

    def analyse(seed):
        rng = np.random.default_rng(seed)
        points = rng.random((2, 100))
        colours = rng.integers(0, 2, size=100, dtype=np.uint16)
        filtration = ch.chromatic.alpha(points, colours, max_num_threads=1)
        return filtration.boundary_matrix()

    with ThreadPoolExecutor(max_workers=4) as pool:
        matrices = list(pool.map(analyse, range(8)))

Here each worker has its own arrays and filtration. ``max_num_threads=1`` keeps
each geometry call serial internally while the workers run separate calls.
Chalc can also use multiple threads within each call. When using a thread pool,
consider limiting the threads within each call to avoid creating more work than
your machine can run at once.

Geometry functions can share input arrays if those arrays stay unchanged. Do not
write to or resize points or colours during a call. This includes changes through
another view of the same NumPy array. Chalc copies the inputs into C++ storage,
but another thread must not change them while that copy is being made.

Sharing a filtration
--------------------

Several threads may read the same ``Filtration`` if no thread changes its
structure, filtration values, or colours. Each thread should create its own
iterator; do not advance the same iterator from several threads.

If the filtration may change, protect both readers and writers with one lock:

.. code-block:: python

    from threading import Lock

    import chalc as ch

    filtration = ch.filtration.complete_complex(4, 2)
    lock = Lock()

    def set_vertex_colour(label, colour):
        with lock:
            vertex = filtration.simplices[0][label]
            vertex.set_colour(colour)
            filtration.propagate_colours()

    def get_boundary_matrix():
        with lock:
            return filtration.boundary_matrix()

The lock covers the colour change and its propagation together, so a reader
cannot see the filtration between those two steps. Locking only the writer would
leave readers unprotected. Changes to an object used by just one thread do not
need a lock.

``Simplex`` objects obtained from a filtration refer to its existing simplices.
The same is true of handles obtained through ``simplices`` or ``facets``.
Reading or changing a simplex must follow the locking rule for its parent
filtration, even if the handle was obtained earlier.

If another thread may add simplices, hold the lock for the whole iteration.
Adding simplices during iteration can invalidate the iterator. Copying the
handles with ``list(filtration)`` does not copy the simplices. Copy the values
you need while holding the lock if you want to use them independently later.

Filtration morphisms
--------------------

``SubChromaticInclusion``, ``KChromaticInclusion``, ``SubChromaticGluing``, and
``KChromaticGluing`` retain the filtration supplied to their constructor.
Creating two morphisms from the same filtration does not make two independent
filtrations.

Concurrent ``sixpack()`` calls may safely share a morphism object if its parameters
and underlying filtration stay unchanged. If either may change, use the same lock
for changing them and for the whole ``sixpack()`` call. Custom subclasses must also
ensure that their own methods and state allow concurrent calls.

Persistence modules and six-packs
---------------------------------

``AugmentedPersistenceModule`` exposes its pairings and entrance times for
reading. Its queries may run concurrently. Do not change its private attributes
or make its entrance-time array writable while other threads use it.

Several threads may read a ``SixPack`` if its data stays unchanged. Its
``SixPackPairings``, ``SixPackRepresentatives``, and ``SixPackMorphisms`` containers
are frozen dataclasses, but their dictionaries and NumPy arrays can still be
changed. ``metadata`` can also be changed, and may contain user objects with
their own threading restrictions. Protect conflicting reads and writes with
the same lock.

``SixPack`` takes ownership of the dataclasses passed to its constructor and
reindexes their contents during construction. Do not use those dataclasses in
another thread while construction is in progress. Passing the same dataclasses
to separate constructors does not create independent data. Make copies first
if you need independent inputs.

Plotting and HDF5 files
-----------------------

Use chalc's plotting and animation functions on the main thread. They use
Matplotlib, whose GUI backends generally require the main thread and whose
objects are not generally safe to access concurrently. See
`Matplotlib's threading guidance
<https://matplotlib.org/stable/users/faq.html#work-with-threads>`_.

``to_hdf5()`` reads the object being saved, so that object must remain unchanged
throughout the call. Separate threads can save to separate files. If threads
share a file or group, use a common lock around each complete save or load,
including opening and closing the file when applicable. A save or load involves
several dataset operations; h5py's internal lock does not make the whole sequence
indivisible. Also prevent changes to a file while another thread is loading it.
See `h5py's threading documentation
<https://docs.h5py.org/en/stable/threads.html>`_.

Free-threaded Python
--------------------

Install chalc and its dependencies for the free-threaded interpreter you use.
An extension that does not declare support may turn the GIL back on when it is
imported. Check after importing chalc and any other libraries your program uses:

.. code-block:: python

    import sys
    import sysconfig

    import chalc

    print(sysconfig.get_config_var("Py_GIL_DISABLED"))  # 1 for a free-threaded build
    print(sys._is_gil_enabled())  # False when the GIL is disabled

The tests in ``tests/test_thread_safety.py`` exercise independent computations,
shared reads, and updates protected by caller locks. They run on regular Python
too. On a free-threaded build they also check that importing chalc leaves the GIL
disabled. These tests cover the supported usage discussed above; they do not make
unsynchronized changes to shared objects safe.
