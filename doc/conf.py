import sphinx_rtd_theme
#import sphinx_fontawesome

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'UNISIS'
copyright = '2024-, Naoto Hori'
author = 'Naoto Hori'
version = 'Draft'
release = 'Draft'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
              'sphinx.ext.mathjax',
              'sphinx.ext.extlinks',
              #'sphinx_fontawesome'
              #'nbsphinx', # for ipynb
]

## MathJax configuration
#mathjax_config = {
#    'TeX': {
#        'equationNumbers': {
#            'autoNumber': 'AMS',
#            'useLabelIds': True
#        }
#    }
#}

#mathjax_options = {
#    'processHtmlClass': 'math',
#    'processEnvironments': True,
#    'processEscapes': True
#}

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
html_static_path = ['_static']

html_css_files = [
    'custom.css',
]

html_js_files = [
    'custom.js',
]

# html_logo  # To show logo

html_theme_options = {
    #'analytics_id': 'G-XXXXXXXXXX',  #  Provided by Google in your dashboard
    #'analytics_anonymize_ip': False,
    #'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': True,
    'vcs_pageview_mode': '',
    #'style_nav_header_background': 'white',

    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
}

html_context = {
    "display_github": True, # Integrate GitHub
    "github_user": "hori-lab", # Username
    "github_repo": "UNISIS", # Repo name
    "github_version": "master", # Version
    "conf_py_path": "/doc/", # Path in the checkout to the docs root
}

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
    "linkify",
        # To automatically identify “bare” web URLs and add hyperlinks:
        #  e.g. www.example.com -> www.example.com
        # To only match URLs that start with schema, such as http://example.com,
        # set myst_linkify_fuzzy_links=False.
    #"replacements",
    #"smartquotes",
    #"strikethrough",
    #"substitution",
    #"tasklist",
]


# -- Options for PDF output -----------------------------------------------
# inside conf.py
#latex_engine = 'xelatex'
latex_elements = {
    'papersize': 'a4paper',
    'pointsize': '11pt',
    'passoptionstopackages': r'''
\PassOptionsToPackage{svgnames}{xcolor}
    ''',
#    'fontpkg': r'''
#\setmainfont{DejaVu Serif}
#\setsansfont{DejaVu Sans}
#\setmonofont{DejaVu Sans Mono}
#    ''',
    'preamble': r'''
\usepackage[titles]{tocloft}
\cftsetpnumwidth {1.25cm}\cftsetrmarg{1.5cm}
\setlength{\cftchapnumwidth}{0.75cm}
\setlength{\cftsecindent}{\cftchapnumwidth}
\setlength{\cftsecnumwidth}{1.25cm}
    ''',
    'sphinxsetup': 'TitleColor=DarkGoldenrod',
    'fncychap': r'\usepackage[Bjornstrup]{fncychap}',
    'printindex': r'\footnotesize\raggedright\printindex',
}
latex_show_urls = 'footnote'
