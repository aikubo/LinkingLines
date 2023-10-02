import os
import sys
sys.path.insert(0, os.path.abspath('../src'))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'LinkingLines'
copyright = '2023, Allison Kubo'
author = 'Allison Kubo'
release = '0.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.mathjax',
              'sphinx.ext.napoleon',
              'sphinx.ext.intersphinx',
              "nbsphinx",
              'myst_parser',
             ]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']
napoleon_google_docstring = False
napoleon_numpy_docstring = True


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
html_theme_options = {
    'github_user': 'aikubo',
    'github_repo': 'Linking-and-Clustering-Dikes',
    'codecov_button': 'true'}
