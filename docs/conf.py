#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# ESSM documentation build configuration file, created by
# sphinx-quickstart on Tue May 16 14:29:58 2017.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

import os
import sys

from sphinx.domains.python import PyModulelevel
from sphinx.ext.autodoc import AutoDirective, ClassDocumenter

from essm.equations._core import BaseEquation
from essm.variables._core import BaseVariable

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
sys.path.insert(0, os.path.abspath('..'))


class BaseDocumenter(ClassDocumenter):
    """Document expression definitions."""

    member_order = 1000

    _BASE_ATTRIBUTES = {
        'default',
        'domain',
        'expr',
        'latex_name',
        'unit',
    }

    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        return isinstance(member, (BaseEquation, BaseVariable)) \
            and hasattr(member, 'definition')

    @property
    def object(self):
        return getattr(self, '_definition', None)

    @object.setter
    def object(self, value):
        if value is not None:
            setattr(self, '_definition', getattr(value, 'definition', value))

    def fiter_members(self, *args, **kwargs):
        members = super(BaseDocumenter, self).filter_members(*args, **kwargs)
        for (mname, member, isattr) in members:
            if mname in self._BASE_ATTRIBUTES:
                yield (mname, member, True)


def process_base(app, objtype, membername, member, skip, options):
    """Always document variables and equations."""
    if isinstance(member, (BaseEquation, BaseVariable)):
        return False


def setup(app):
    app.connect('autodoc-skip-member', process_base)
    AutoDirective._registry['base_expression'] = BaseDocumenter

    def get_attr(obj, value, *args, **kwargs):
        return getattr(obj.definition, value)

    AutoDirective._special_attrgetters[BaseVariable] = get_attr
    AutoDirective._special_attrgetters[BaseEquation] = get_attr


# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Do not warn on external images.
suppress_warnings = ['image.nonlocal_uri', 'ref.citation']

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'matplotlib.sphinxext.plot_directive', 'sphinxcontrib.bibtex',
    'sphinx.ext.autodoc', 'sphinx.ext.autosummary', 'sphinx.ext.doctest',
    'sphinx.ext.intersphinx', 'sphinx.ext.todo', 'sphinx.ext.coverage',
    'sphinx.ext.mathjax', 'sphinx.ext.viewcode', 'sphinx.ext.githubpages',
    'nbsphinx'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = 'ESSM'
copyright = '2017, Stan Schymanski'
author = 'Stan Schymanski'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.

# Get the version string. Cannot be done with import!
g = {}
with open(os.path.join(os.path.dirname(__file__), '..', 'essm', 'version.py'),
          'rt') as fp:
    exec(fp.read(), g)
    version = g['__version__']

# The full version, including alpha/beta/rc tags.
release = version

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store',
                   '**.ipynb_checkpoints']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
    'logo': 'logo.png',
    'logo_name': False,
    'description':
        'Python package for dealing with physical variables '
        'and units.',
    'github_user': 'environmentalscience',
    'github_repo': 'essm',
    'github_button': False,
    'github_banner': True,
    'show_powered_by': False,
    'extra_nav_links': {
        'essm@GitHub': 'https://github.com/environmentalscience/essm',
        'essm@PyPI': 'https://pypi.python.org/pypi/essm/',
    }
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom sidebar templates, maps document names to template names.
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',
        'searchbox.html',
        'donate.html',
    ]
}

# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'ESSMdoc'

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (
        master_doc, 'ESSM.tex', 'ESSM Documentation', 'Stan Schymanski',
        'manual'
    ),
]

# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, 'essm', 'ESSM Documentation', [author], 1)]

# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc, 'ESSM', 'ESSM Documentation', author, 'ESSM',
        'One line description of project.', 'Miscellaneous'
    ),
]

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'https://docs.python.org/': None}

autoclass_content = 'both'

autodoc_default_options = {
    'undoc-members': True,
    'show-inheritance': True
}
