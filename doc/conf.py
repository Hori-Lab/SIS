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

extensions = [
              'sphinx.ext.mathjax',
              #'sphinx_fontawesome'
              #'nbsphinx', # for ipynb
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '_sphinx.md', '.git']

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
    #'.txt': 'markdown',
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
#html_theme = 'alabaster'

extensions += ['sphinx_rtd_theme',]
html_theme = 'sphinx_rtd_theme'
#html_static_path = ['_static']

# -- MyST --------------------------------------------------------------------
extensions += ['myst_parser',]
myst_enable_extensions = [
    "amsmath",
    #"colon_fence",
    #"deflist",
    "dollarmath",
    #"fieldlist",
    #"html_admonition",
    #"html_image",
    #"linkify",
    #"replacements",
    #"smartquotes",
    #"strikethrough",
    #"substitution",
    #"tasklist",
]

