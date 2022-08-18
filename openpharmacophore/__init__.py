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
from .pharmacophore.pharmacophoric_point import PharmacophoricPoint
from .pharmacophore.pharmacophore import Pharmacophore
from .pharmacophore.ligand_based import LigandBasedPharmacophore
from .pharmacophore.structured_based import StructuredBasedPharmacophore
from .screening.screening import VirtualScreening
from .screening.multiprocess import MultiProcessVirtualScreening
from .screening.retrospective import RetrospectiveScreening
from .pharmacophore.dynophore import Dynophore
from .databases.zinc_client import ZincClient
