import os
import sys
sys.path.insert(0, os.path.abspath('../../'))
sys.path.insert(0, os.path.abspath('../../src/'))

sys.path.insert(0, os.path.abspath('../../src/linkinglines/'))
sys.path.insert(0, os.path.abspath('../../src/linkinglines/HT'))
sys.path.insert(0, os.path.abspath('../../src/linkinglines/ProcessingUtils'))
sys.path.insert(0, os.path.abspath('../../src/linkinglines/FeatureExtraction'))
sys.path.insert(0, os.path.abspath('../../src/linkinglines/PlotUtils'))
sys.path.insert(0, os.path.abspath('../../src/linkinglines/ClusterLines'))


import linkinglines


# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'LinkingLines'
copyright = '2023, Allison Kubo'
author = 'Allison Kubo'
release = '2.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.mathjax',
              'sphinx.ext.napoleon',
              'sphinx.ext.intersphinx',
              'myst_parser',
             ]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# The encoding of source files.
# source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = "index"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
html_theme_options = {
    'github_user': 'aikubo',
    'github_repo': 'Linking-and-Clustering-Dikes',
    "github_banner": True,
    "show_related": False,
    "note_bg": "#FFF59C"}
