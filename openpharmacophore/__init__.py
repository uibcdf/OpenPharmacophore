"""
OpenPharmacophore
A library to work with pharmacophores.
"""

# Handle versioneer
from ._version import get_versions


import pyunitwizard as puw

puw.configure.load_library(['pint', "openmm.unit"])
puw.configure.set_default_form('pint')
puw.configure.set_standard_units(['angstroms', 'ps', 'K', 'mole', 'amu', 'e',
                                 'kJ/mol', 'kJ/(mol*nm**2)', 'N', 'degrees'])


# Add imports here
from .molecular_systems import Ligand, LigandSet, Protein, BindingSite
from .point.pharmacophoric_point import PharmacophoricPoint, distance_between_pharmacophoric_points
from .pharmacophore.pharmacophore import Pharmacophore
from .pl_complex.pl_complex import PLComplex
from .pharmacophore.ligand_based.ligand_based import LigandBasedPharmacophore
from .pharmacophore.ligand_receptor.ligand_receptor import LigandReceptorPharmacophore
from .load.load import load, load_ligands
from .screening import VirtualScreening


versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

__documentation_web__ = 'https://www.uibcdf.org/openpharmacophore'
__github_web__ = 'https://github.com/uibcdf/openpharmacophore'
__github_issues_web__ = __github_web__ + '/issues'
