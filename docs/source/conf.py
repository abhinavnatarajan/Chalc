import sys
from pathlib import Path

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here.
stubs_path = (Path(__file__).resolve().parents[2] / 'stubs' ).resolve()
sys.path.insert(0, str(stubs_path))
sys.path.append(str((Path(__file__).resolve().parents[1] / 'exts').resolve()))

project = 'Chalc'
copyright = '2023, Abhinav Natarajan'
author = 'Abhinav Natarajan'
release = '0.3.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

maximum_signature_line_length = 80
toc_object_entries_show_parents = 'hide'
add_module_names = False

extensions = [
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.napoleon',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'autosummary_customise',
    'sphinx.ext.intersphinx',
    'sphinx_toolbox.sidebar_links',
    'sphinx_toolbox.github',
    'sphinx.ext.githubpages']

# autosectionlabel options
autosectionlabel_prefix_document = True

# napolean options
# napoleon_include_init_with_doc = True
# napoleon_preprocess_types = True

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

master_doc = "index"
templates_path = ['templates']
# exclude_patterns = ['']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_static_path = ['_static']
html_css_files = [
    'css/custom.css',
]

html_theme = 'furo'
html_theme_options = {
    "source_repository": "https://github.com/abhinavnatarajan/chalc",
    "source_branch": "master",
    "source_directory": "docs/",
    "top_of_page_button": "edit"
}