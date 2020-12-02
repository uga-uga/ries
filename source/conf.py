import os
import sys
sys.path.insert(0, os.path.abspath('../ries/'))
sys.path.insert(0, os.path.abspath('../ries/constituents'))
sys.path.insert(0, os.path.abspath('../ries/integration'))
sys.path.insert(0, os.path.abspath('../ries/nonresonant'))


# -- Project information -----------------------------------------------------

project = 'ries'
copyright = '2020, Udo Friman-Gayer'
author = 'Udo Friman-Gayer'


# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinxcontrib.bibtex'
]

# -- Options for HTML output -------------------------------------------------

html_theme = 'classic'