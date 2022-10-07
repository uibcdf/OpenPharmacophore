"""
OpenPharmacophore
A library to work with pharmacophores.
"""

# Handle versioneer
from ._version import get_versions

# Add imports here
from ._pyunitwizard import puw
from .pharmacophore import PharmacophoricPoint, distance_between_pharmacophoric_points
from .pharmacophore import LigandBasedPharmacophore
from .pharmacophore import StructureBasedPharmacophore
from .screening import VirtualScreening



versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

__documentation_web__ = 'https://www.uibcdf.org/openpharmacophore'
__github_web__ = 'https://github.com/uibcdf/openpharmacophore'
__github_issues_web__ = __github_web__ + '/issues'
