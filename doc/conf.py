import os

from setuptools_scm import get_version

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
# sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, "src")))

extensions = [
    "sphinx.ext.autodoc",
    "sphinx_issues",
    "sphinx_better_subsection",
]

source_suffix = ".rst"
master_doc = "index"
project = "Cutadapt"
copyright = "2010 Marcel Martin"

version = get_version(root="..", relative_to=__file__)

# Read The Docs modifies the conf.py script and we therefore get
# version numbers like 0.12+0.g27d0d31
if os.environ.get("READTHEDOCS") == "True":
    version = ".".join(version.split(".")[:2])
    html_theme = "sphinx_rtd_theme"

# The full version, including alpha/beta/rc tags.
release = version

issues_uri = "https://github.com/marcelm/cutadapt/issues/{issue}"
issues_pr_uri = "https://github.com/marcelm/cutadapt/pull/{pr}"
suppress_warnings = ["image.nonlocal_uri"]
exclude_patterns = ["_build"]
pygments_style = "sphinx"

html_static_path = ["_static"]
