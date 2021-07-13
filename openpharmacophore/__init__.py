"""
OpenPharmacophore
This must be a short description of the project
"""

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

__documentation_web__ = 'https://www.uibcdf.org/OpenPharmacophore'
__github_web__ = 'https://github.com/uibcdf/OpenPharmacophore'
__github_issues_web__ = __github_web__ + '/issues'

# Add imports here
from ._pyunitwizard import puw as _puw
from . import pharmacophoric_elements
from .pharmacophore import Pharmacophore
from . import demo as demo
