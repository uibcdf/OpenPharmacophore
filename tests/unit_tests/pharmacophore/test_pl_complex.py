import openpharmacophore.data as data
from openpharmacophore.pharmacophore.pl_complex import PLComplex
import pytest


@pytest.fixture()
def pl_complex():
    """ Returns a pl complex for testing. """
    return PLComplex(data.pdb["binding_site.pdb"])


def test_init_pl_complex(pl_complex):
    assert pl_complex.traj.n_atoms == 166
    assert pl_complex.topology.n_chains == 2
    assert pl_complex.mol_graph.GetNumAtoms() == 166


def test_is_ligand_atom_protein_atom(mocker):
    mock_atom = mocker.Mock()
    mock_atom.residue.is_water = False
    mock_atom.residue.is_protein = True
    mock_atom.residue.n_atoms = 8
    assert not PLComplex._is_ligand_atom(mock_atom)


def test_is_ligand_atom_water_atom(mocker):
    mock_atom = mocker.Mock()
    mock_atom.residue.is_water = True
    mock_atom.residue.is_protein = False
    mock_atom.residue.n_atoms = 1
    assert not PLComplex._is_ligand_atom(mock_atom)


def test_is_ligand_atom_ligand_atom(mocker):
    mock_atom = mocker.Mock()
    mock_atom.residue.is_water = False
    mock_atom.residue.is_protein = False
    mock_atom.residue.n_atoms = 15
    assert PLComplex._is_ligand_atom(mock_atom)


def test_find_ligand(pl_complex):
    assert pl_complex.find_ligands() == ["EST:B"]
