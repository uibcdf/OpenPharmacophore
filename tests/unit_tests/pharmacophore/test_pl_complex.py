import openpharmacophore._private_tools.exceptions as exc
import openpharmacophore.data as data
from openpharmacophore.pharmacophore.pl_complex import PLComplex
import pytest


@pytest.fixture()
def pl_complex():
    """ Returns a pl complex for testing. """
    return PLComplex(data.pdb["test_with_lig.pdb"])


def test_init_pl_complex(pl_complex):
    assert pl_complex.traj.n_atoms == 166
    assert pl_complex.topology.n_chains == 2
    assert pl_complex.mol_graph.GetNumAtoms() == 166


def test_init_pl_complex_with_no_ligand_raises_error():
    with pytest.raises(exc.NoLigandError):
        PLComplex(data.pdb["test_no_lig.pdb"])


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


def test_find_ligands(pl_complex):
    assert pl_complex.ligand_ids == ["EST:B"]


def test_ligand_indices(pl_complex):
    pl_complex._ligand_atom_indices("EST:B")
    expected_indices = list(range(146, 166))
    assert len(pl_complex._lig_indices) == len(expected_indices)
    assert pl_complex._lig_indices == expected_indices


def test_ligand_to_mol(pl_complex):
    pl_complex._ligand_atom_indices("EST:B")
    pl_complex._ligand_to_mol()
    assert pl_complex.ligand.GetNumAtoms() == 20


def test_ligand_to_mol_empty_indices_list_raises_error(pl_complex):
    with pytest.raises(exc.NoLigandIndicesError):
        pl_complex._ligand_to_mol()
