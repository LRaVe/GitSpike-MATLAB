import os
import sys

project = 'Discriminative neuronal subpopulation algorithms'
copyright = '2026, Maxime BELTOISE & Laure WOLFF'
author = 'Maxime BELTOISE & Laure WOLFF'

# Chemin direct vers Summed_Population
summed_pop_dir = os.path.abspath('../../Summed_Population')
matlab_src_dir = summed_pop_dir

# On ajoute le dossier et ses sous-dossiers au path
sys.path.insert(0, summed_pop_dir)
for dirpath, dirnames, filenames in os.walk(summed_pop_dir):
    sys.path.append(dirpath)

primary_domain = 'mat'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinxcontrib.matlab',
]

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']