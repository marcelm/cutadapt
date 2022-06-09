import sys
import os
import time
from pkg_resources import get_distribution

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, "src")))

extensions = [
    "sphinx.ext.autodoc",
    "sphinx_issues",
    "sphinx_better_subsection",
]

source_suffix = ".rst"
master_doc = "index"
project = "Cutadapt"
copyright = "2010-{}, Marcel Martin".format(time.gmtime().tm_year)

release = get_distribution("cutadapt").version
# Read The Docs modifies the conf.py script and we therefore get
# version numbers like 0.12+0.g27d0d31
if os.environ.get("READTHEDOCS") == "True":
    version = ".".join(release.split(".")[:2])
else:
    version = release

# The full version, including alpha/beta/rc tags.
release = version

issues_uri = "https://github.com/marcelm/cutadapt/issues/{issue}"
issues_pr_uri = "https://github.com/marcelm/cutadapt/pull/{pr}"
suppress_warnings = ["image.nonlocal_uri"]
exclude_patterns = ["_build"]
pygments_style = "sphinx"

try:
    from better import better_theme_path

    html_theme_path = [better_theme_path]
    html_theme = "better"
except ImportError:
    html_theme = "default"

html_static_path = ["_static"]
