project = 'ogtk'
copyright = '2025, Pedro Olivares'
author = 'Pedro Olivares'

extensions = [
    # moving away from mk for docs
    # https://www.ericholscher.com/blog/2016/mar/15/dont-use-markdown-for-technical-docs/
    #'myst_parser',
    'sphinx.ext.duration',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',  # Optional: for Google/NumPy-style docstrings
    'sphinx.ext.viewcode',  # Optional: for linking to source code
]

source_suffix = {
    '.rst': 'restructuredtext',
    #'.md': 'markdown',
        }

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

#html_theme = 'alabaster'
html_theme = 'furo'
#html_theme = 'sphinx_book_theme'
html_static_path = ['_static']

#html_logo = "folded.jpg"
html_logo = "honeyomics4.jpg"

autosummary_generate = True

from sphinx_pyproject import SphinxConfig
config = SphinxConfig("../pyproject.toml", globalns=globals())
