# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
from generate_logo import make_logo
import os
import sys
import datetime
import re
import subprocess
sys.path.insert(0, os.path.abspath('.'))

DOCS_SOURCE_DIR = os.path.abspath(os.path.dirname(__file__))
REPO_ROOT_DIR = os.path.abspath(os.path.join(DOCS_SOURCE_DIR, '..', '..'))
BUILD_DIR = os.path.join(REPO_ROOT_DIR, 'build')
DOXYGEN_XML_DIR = os.path.join(BUILD_DIR, 'docs', 'doxygen', 'xml')


def _run_command(command, cwd=None, env=None):
    try:
        subprocess.run(command, cwd=cwd, env=env, check=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        sys.stderr.write(f"Command failed: {command}\n{exc}\n")
        raise


def _build_cmake_env():
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if not conda_prefix:
        return None

    env = os.environ.copy()
    pkg_config_path = os.path.join(conda_prefix, 'lib', 'pkgconfig')
    existing_pkg_config_path = env.get('PKG_CONFIG_PATH')
    if existing_pkg_config_path:
        env['PKG_CONFIG_PATH'] = f"{pkg_config_path}:{existing_pkg_config_path}"
    else:
        env['PKG_CONFIG_PATH'] = pkg_config_path

    existing_cmake_prefix_path = env.get('CMAKE_PREFIX_PATH')
    if existing_cmake_prefix_path:
        env['CMAKE_PREFIX_PATH'] = f"{conda_prefix}:{existing_cmake_prefix_path}"
    else:
        env['CMAKE_PREFIX_PATH'] = conda_prefix

    return env


def _ensure_doxygen_xml():
    read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'
    if not read_the_docs_build:
        return

    cmake_env = _build_cmake_env()

    cache_file = os.path.join(BUILD_DIR, 'CMakeCache.txt')
    if not os.path.exists(cache_file):
        _run_command([
            'cmake',
            '-S', REPO_ROOT_DIR,
            '-B', BUILD_DIR,
            '-D', 'TUVX_DOCS_ONLY=ON'
        ], env=cmake_env)

    _run_command([
        'cmake',
        '--build', BUILD_DIR,
        '--target', 'Doxygen'
    ], env=cmake_env)


# -- Project information -----------------------------------------------------
project = 'TUV-x'
copyright = f"2022-{datetime.datetime.now().year}, NSF-NCAR/ACOM"
author = 'NSF-NCAR/ACOM'

# the suffix is required. This is controlled by the dockerfile that builds the docs
regex = r'project\(.*VERSION\s+(\d+\.\d+\.\d+)'
version = '0.0.0'
# read the version from the cmake files
with open(f'../../CMakeLists.txt', 'r') as f:
    for line in f:
        match = re.match(regex, line)
        if match:
            version = match.group(1)
release = f'v{version}'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_design',
    'sphinx.ext.todo',
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'restructuredtext',
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

todo_include_todos = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'pydata_sphinx_theme'

html_theme_options = {
    "external_links": [],
    "github_url": "https://github.com/NCAR/tuv-x",
    "navbar_end": ["navbar-icon-links"],
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

make_logo(tuvx_static_path=html_static_path[0])

# https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-html_css_files
# Custom styling
html_css_files = [
    'custom.css'
]

html_favicon = '_static/favicon.ico'

html_logo = '_static/logo.svg'


def setup(app):
    app.connect("builder-inited", lambda _app: _ensure_doxygen_xml())
