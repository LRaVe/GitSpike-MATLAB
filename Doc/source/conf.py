# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

_HERE = os.path.dirname(__file__)
_ROOT = os.path.abspath(os.path.join(_HERE, '..', '..'))

sys.path.insert(0, _ROOT)

project = 'GitSpike-MATLAB'
copyright = '2026, Lucas RAVELOARINORO, Agathe JULIEN, Laure WOLFF, Maxime BELTOISE'
author = 'Lucas RAVELOARINORO, Agathe JULIEN, Laure WOLFF, Maxime BELTOISE'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc','sphinxcontrib.matlab','sphinx.ext.viewcode','sphinx.ext.napoleon']
primary_domain = 'mat'
matlab_src_dir = _ROOT

templates_path = ['_templates']
exclude_patterns = ['_generated']

html_sidebars = {
	'**': ['localtoc.html', 'relations.html', 'searchbox.html'],
}


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
