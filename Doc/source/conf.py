# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# -- Path setup --------------------------------------------------------------

_HERE = os.path.dirname(__file__)
_ROOT = os.path.abspath(os.path.join(_HERE, '..', '..'))

# Ajout prioritaire de la racine et du dossier ALGO_Simple au path Python
sys.path.insert(0, _ROOT)
algo_simple_dir = os.path.join(_ROOT, 'ALGO_Simple')
sys.path.insert(0, algo_simple_dir)

# Parcourt récursivement tous les sous-dossiers du projet
for root, dirs, files in os.walk(_ROOT):
    for d in dirs:
        sys.path.insert(0, os.path.join(root, d))


# -- Project information -----------------------------------------------------

project = 'GitSpike-MATLAB'
copyright = '2026, Lucas RAVELOARINORO, Agathe JULIEN, Laure WOLFF, Maxime BELTOISE'
author = 'Lucas RAVELOARINORO, Agathe JULIEN, Laure WOLFF, Maxime BELTOISE'
release = '0.1.0'


# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinxcontrib.matlab',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon'
]

primary_domain = 'mat'
matlab_src_dir = _ROOT

# Configuration de Napoleon pour extraire correctement les blocs 'Parameters' et 'Returns'
napoleon_custom_sections = [('Returns', 'params_style')]

templates_path = ['_templates']
exclude_patterns = ['_generated', 'Thumbs.db', '.DS_Store']

matlab_auto_link = False


# -- Correctif pour le bug 'NoneType' object has no attribute 'args' --
def skip_bad_signatures(app, what, name, obj, options, signature, return_annotation):
    # Intercepte l'erreur de signature pour éviter que Sphinx 7.4 ne plante sur le formatage
    return signature, return_annotation

def setup(app):
    # Désactive l'inspection automatique problématique des arguments
    app.connect('autodoc-process-signature', skip_bad_signatures)


# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

html_sidebars = {
    '**': ['localtoc.html', 'relations.html', 'searchbox.html'],
}