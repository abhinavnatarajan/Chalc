import os
from pathlib import Path

from packaging.version import parse

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here.
docs_path = list(Path(__file__).parents)[1]
project_root_dir = docs_path.parent

copyright = "2023, Abhinav Natarajan"
author = "Abhinav Natarajan"

# Get the name and release from pyproject.toml
from tomllib import loads

with open(project_root_dir / "pyproject.toml") as f:
	proj_props = loads(f.read())["project"]

project = proj_props["name"]
release = ".".join(map(str, parse(proj_props["version"]).release))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

maximum_signature_line_length = 160
toc_object_entries_show_parents = "hide"
add_module_names = False
master_doc = "index"
# templates_path = ['templates']
exclude_patterns = [
	"_build",
	"templates",
	"exts",
	"static",
	# "example.ipynb",
]

extensions = [
	"sphinx.ext.autosectionlabel",
	"sphinx.ext.napoleon",
	"autoapi.extension",
	"sphinx.ext.intersphinx",
	"nbsphinx",
	"sphinx_design",
	"sphinx_toolbox.sidebar_links",
	"sphinx_toolbox.github",
	"sphinx.ext.githubpages",
]

# autosectionlabel options
autosectionlabel_prefix_document = True

# napoleon options
napolean_include_init_with_doc = True
napoleon_attr_annotations = True

# autoapi options
autoapi_dirs = ["../../src"]
autoapi_options = [
	"members",
	"show-module-summary",
	"special-members",
	"inherited-members",
	"show-inheritance",
]
autoapi_member_order = "groupwise"

# intersphinx options
intersphinx_mapping = {
	"python": ("https://docs.python.org/3/", None),
	"numpy": ("https://numpy.org/doc/stable", None),
	"matplotlib": ("https://matplotlib.org/stable", None),
	"h5py": ("https://docs.h5py.org/en/stable/", None),
}

# nbsphinx options
nbsphinx_execute = "always"
os.environ["PYDEVD_DISABLE_FILE_VALIDATION"] = "1"

# toolbox.github options
github_username = "abhinavnatarajan"
github_repository = "chalc"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_static_path = ["../static"]
html_css_files = ["css/custom.css"]

html_theme = "furo"
html_theme_options = {
	"source_repository": "https://github.com/abhinavnatarajan/chalc",
	"source_branch": "master",
	"source_directory": "docs/",
	"top_of_page_button": "edit",
}
