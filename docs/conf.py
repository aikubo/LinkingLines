import os
import sys
sys.path.insert(0, os.path.abspath('../src/linkinglines'))

sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.abspath('./notebooks'))


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
release = '0.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.mathjax',
              'sphinx.ext.napoleon',
              'sphinx.ext.intersphinx',
              'nbsphinx']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
napoleon_google_docstring = False
napoleon_numpy_docstring = True
source_suffix = [".rst", ".md"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

# Add GitHub badge to the sidebar
html_theme_options = {'description': "Using the Hough Transform to cluster lines",
                      'github_user': 'aikubo',
                      'github_repo': 'LinkingLines',
                      'github_banner': True,
}


html_logo = "_static/dikeslogo.svg"
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',
        'searchbox.html',
        
    ]
}

