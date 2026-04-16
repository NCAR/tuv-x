import datetime
import os
import re
import subprocess
import sys

DOCS_SOURCE_DIR = os.path.abspath(os.path.dirname(__file__))
REPO_ROOT_DIR = os.path.abspath(os.path.join(DOCS_SOURCE_DIR, '..', '..'))
BUILD_DIR = os.path.join(REPO_ROOT_DIR, 'build')
DOXYGEN_XML_DIR = os.path.join(BUILD_DIR, 'docs', 'doxygen', 'xml')


def _run_command(command, cwd=None):
    try:
        subprocess.run(command, cwd=cwd, check=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        sys.stderr.write(f"Command failed: {command}\n{exc}\n")
        raise


def _ensure_doxygen_xml():
    # On Read the Docs, drive the CMake Doxygen target so breathe has XML to consume.
    read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'
    if not read_the_docs_build:
        return

    cache_file = os.path.join(BUILD_DIR, 'CMakeCache.txt')
    if not os.path.exists(cache_file):
        _run_command([
            'cmake',
            '-S', REPO_ROOT_DIR,
            '-B', BUILD_DIR,
            '-D', 'TUVX_BUILD_DOCS=ON',
            '-D', 'TUVX_ENABLE_TESTS=OFF',
        ])

    _run_command([
        'cmake',
        '--build', BUILD_DIR,
        '--target', 'Doxygen',
    ])


project = 'TUV-x'
copyright = f'2020-{datetime.datetime.now().year}, UCAR'
author = 'NSF-NCAR/ACOM'

version = '0.0.0'
with open(os.path.join(REPO_ROOT_DIR, 'CMakeLists.txt'), 'r') as f:
    for line in f:
        match = re.match(r'project\([\w-]+\s+VERSION\s+(\d+\.\d+\.\d+)', line)
        if match:
            version = match.group(1)
            break
release = version

extensions = [
    'breathe',
    'sphinx_copybutton',
    'sphinx_design',
]

breathe_default_project = 'tuvx'
breathe_projects = {
    'tuvx': DOXYGEN_XML_DIR,
}

exclude_patterns = []

html_theme = 'pydata_sphinx_theme'
html_theme_options = {
    'github_url': 'https://github.com/NCAR/tuv-x',
    'navbar_end': ['navbar-icon-links'],
}


def setup(app):
    app.connect('builder-inited', lambda _app: _ensure_doxygen_xml())
