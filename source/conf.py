# This file is part of ries.

# ries is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ries is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ries.  If not, see <https://www.gnu.org/licenses/>.

import os
import sys
sys.path.insert(0, os.path.abspath('../ries/'))
sys.path.insert(0, os.path.abspath('../ries/constituents'))
sys.path.insert(0, os.path.abspath('../ries/integration'))
sys.path.insert(0, os.path.abspath('../ries/nonresonant'))
sys.path.insert(0, os.path.abspath('../ries/resonance'))

# -- Project information -----------------------------------------------------

project = 'ries'
copyright = '2020, Udo Friman-Gayer'
author = 'Udo Friman-Gayer'


# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinxcontrib.bibtex'
]
bibtex_bibfiles = ['bibliography.bib']

# -- Options for HTML output -------------------------------------------------

html_theme = 'classic'
