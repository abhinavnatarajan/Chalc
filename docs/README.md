# Chalc Documentation Source Directory

To build the documentation, first make sure that ``chalc`` is installed and available to the Python interpreter. Then run the following commands from the root directory of the repository.

```bash
pip install -r requirements.txt
make html
```

This will build the documentation into the folder ``docs/build`` with root ``index.html``.
In order to avoid polluting your global Python package tree, it is a good idea to install ``chalc`` and run the ``pip install ..`` command above in a virtual environment.
