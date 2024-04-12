import sphinx_rtd_theme
#import sphinx_fontawesome

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SIS'
copyright = '2024, Naoto Hori'
author = 'Naoto Hori'
release = '2024.04'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['myst_parser',
              'sphinx_rtd_theme',
              #'sphinx_fontawesome'
              'sphinx.ext.mathjax'
            ]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '_sphinx.md', '.git']

source_suffix = {
    '.rst': 'restructuredtext',
    #'.txt': 'markdown',
    '.md': 'markdown',
}

mathjax_path = '.'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
