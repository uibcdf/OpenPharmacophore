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
from .molecular_systems import Ligand, Protein, smiles_from_pdb_id
from .molecular_systems import ComplexBindingSite, BindingSite
from .pharmacophore.pharmacophoric_point import PharmacophoricPoint
from .pharmacophore.pharmacophore import Pharmacophore
from .pharmacophore.ligand_based.ligand_based import LigandBasedPharmacophore
from .pharmacophore.ligand_receptor.ligand_receptor import LigandReceptorPharmacophore
from .visualization import Viewer, draw_ligands, draw_ligands_chem_feats
from .io import mol_files, pharmacophore_reader, pharmacophore_writer
from .load import load, load_ligands
from .screening.mol_db import MolDB


versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

__documentation_web__ = 'https://www.uibcdf.org/openpharmacophore'
__github_web__ = 'https://github.com/uibcdf/openpharmacophore'
__github_issues_web__ = __github_web__ + '/issues'
