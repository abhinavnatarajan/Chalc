import os
from pathlib import Path
from tomllib import loads

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

copyright = "2023, Abhinav Natarajan"  # noqa: A001
author = "Abhinav Natarajan"

# Get the name and release from pyproject.toml

with Path.open(project_root_dir / "pyproject.toml") as f:
	proj_props = loads(f.read())["project"]

project = proj_props["name"]
release = ".".join(map(str, parse(proj_props["version"]).release))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

maximum_signature_line_length = 160
toc_object_entries_show_parents = "hide"
add_module_names = False
master_doc = "index"
exclude_patterns = [
	"build",
	"static",
	# "example.ipynb",
]

extensions = [
	"autoapi.extension",
	"sphinx.ext.autosectionlabel",
	"sphinx.ext.napoleon",
	"sphinx.ext.intersphinx",
	"nbsphinx",
	"sphinx_design",
	"sphinx_toolbox.sidebar_links",
	"sphinx_toolbox.github",
	"sphinx.ext.githubpages",
	"sphinx.ext.inheritance_diagram",
]

# autosectionlabel options
autosectionlabel_prefix_document = True

# napoleon options
napolean_include_init_with_doc = True
napoleon_attr_annotations = True

# autoapi options
autoapi_dirs = ["../../src"]
autoapi_add_toctree_entry = False
autoapi_options = [
	"members",
	"show-module-summary",
	"special-members",
	# "inherited-members",
	"show-inheritance",
	"show-inheritance-diagram",
	"imported-members",
]
autoapi_member_order = "alphabetical"
autoapi_keep_files = True
def autoapi_skip_member(app, what, name, obj, skip, options):
	if name == "chalc.sixpack.types.DiagramName":
		skip = True
	return skip

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

def setup(sphinx):
	sphinx.connect("autoapi-skip-member", autoapi_skip_member)
