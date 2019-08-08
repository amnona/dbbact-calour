# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# import sphinx.ext.autosummary as autosummary
import sphinx_bootstrap_theme

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# ignore lots of module imports to make readthedocs environment
autodoc_mock_imports = ['numpy', 'scipy', 'pandas', 'calour', 'statsmodels', 'sip', 'PyQt5', 'PyQt5.QtGui', 'PyQt5.QtCore', 'PyQt5.QtWidgets']

extensions = [
    'sphinx.ext.autodoc',
    # 'sphinx_autodoc_typehints',   # something wrong with the latest version 1.2.5
    'sphinx.ext.mathjax',
    # 'numpydoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    # 'nbsphinx',
    # this extenstion is needed to avoid the
    # "WARNING: Pygments lexer name 'ipython3' is not known" error for the notebooks.
    # fix based on: https://github.com/spatialaudio/nbsphinx/issues/24
    # 'IPython.sphinxext.ipython_console_highlighting'
]


# The suffix of source filenames.
source_suffix = '.rst'

# The encoding of source files.
# source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'


# -- Options for napoleon -------------------------------------------------

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = False
napoleon_use_ivar = True
napoleon_use_param = True
napoleon_use_rtype = True


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'bootstrap'
html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()

# Theme options are theme-specific and customize the look and feel of a theme
# further. For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    # disable the sidebar
    'nosidebar': True,

    # Navigation bar title. (Default: ``project`` value)
    'navbar_title': 'dbbact-calour docs',

    # Render the next and previous page links in navbar. (Default: true)
    'navbar_sidebarrel': False,

    # Bootswatch (http://bootswatch.com/) theme.
    #
    # Options are nothing with "" (default) or the name of a valid theme
    # such as "amelia" or "cosmo".
    'bootswatch_theme': 'united',

    # Location of link to source.
    # Options are "nav" (default), "footer" or anything else to exclude.
    'source_link_position': False,

    # Choose Bootstrap version.
    'bootstrap_version': "3"
}

# -- Options for autosummary ----------------------------------------------
autosummary_generate = True

# -- Project information -----------------------------------------------------

project = 'dbbact-calour'
copyright = '2019, dbbact calour team'
author = 'dbbact calour team'


# -- General configuration ---------------------------------------------------

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# # -- Options for HTML output -------------------------------------------------

# # The theme to use for HTML and HTML Help pages.  See the documentation for
# # a list of builtin themes.
# #
# html_theme = 'alabaster'

# # Add any paths that contain custom static files (such as style sheets) here,
# # relative to this directory. They are copied after the builtin static files,
# # so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']
