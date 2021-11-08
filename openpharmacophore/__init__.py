"""
OpenPharmacophore
A library to work with pharmacophores.
"""

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

__documentation_web__ = 'https://www.uibcdf.org/openpharmacophore'
__github_web__ = 'https://github.com/uibcdf/openpharmacophore'
__github_issues_web__ = __github_web__ + '/issues'

# Add imports here
from ._pyunitwizard import puw as _puw
from .pharmacophoric_point import PharmacophoricPoint
from .pharmacophore import Pharmacophore
from .ligand_based import LigandBasedPharmacophore
from .structured_based import StructuredBasedPharmacophore
from .screening.screening import VirtualScreening
from .screening.retrospective import RetrospectiveScreening
from .dynophore import Dynophore as Dynophore
from . import demo as demo
