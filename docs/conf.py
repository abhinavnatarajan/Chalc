import sys
from pathlib import Path
from platform import python_version
from packaging.version import parse

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here.
docs_path = Path(__file__).parent
project_root_dir = docs_path.parent
stubs_path = str((docs_path / 'stubs' ).resolve())
exts_path = str((docs_path / 'exts').resolve())
sys.path.insert(0, stubs_path)
sys.path.append(exts_path)

copyright = '2023, Abhinav Natarajan'
author = 'Abhinav Natarajan'

# Get the name and release from pyproject.toml
if parse(python_version()) >= parse('3.11'):
    from tomllib import loads as tomlread
else:
    from toml import loads as tomlread
with open(project_root_dir / "pyproject.toml") as f:
    proj_props = tomlread(f.read())['project']

project = proj_props['name']
release = '.'.join(map(str, parse(proj_props['version']).release))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

maximum_signature_line_length = 80
toc_object_entries_show_parents = 'hide'
add_module_names = False
master_doc = "index"
templates_path = ['_templates']
# exclude_patterns = ['example.ipynb']

extensions = [
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.napoleon',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'autosummary_customise',
    'sphinx.ext.intersphinx',
    'nbsphinx',
    'sphinx_design',
    'sphinx_toolbox.sidebar_links',
    'sphinx_toolbox.github',
    'sphinx.ext.githubpages']

# autosectionlabel options
autosectionlabel_prefix_document = True

# autodoc options
autodoc_typehints = 'description'
autodoc_typehints_description_target = "documented_params"
autodoc_mock_imports = ['numpy']

#autosummary options
autosummary_generate = True

# intersphinx options
intersphinx_mapping = { 
    'python' : ('https://docs.python.org/3/', None),
    'numpy' : ('https://numpy.org/doc/stable', None)}

# toolbox.github options
github_username = 'abhinavnatarajan'
github_repository = 'chalc'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_static_path = ['_static']
html_css_files = [
    'css/custom.css'
]

html_theme = 'furo'
html_theme_options = {
    "source_repository": "https://github.com/abhinavnatarajan/chalc",
    "source_branch": "master",
    "source_directory": "docs/",
    "top_of_page_button": "edit"
}