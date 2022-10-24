import openpharmacophore._private_tools.exceptions as exc
import openpharmacophore.data as data
from openpharmacophore.pharmacophore.pl_complex import PLComplex
import pytest
from copy import deepcopy
from rdkit import Chem


@pytest.fixture()
def pl_complex():
    """ Returns a pl complex for testing. """
    pl_complex = PLComplex(data.pdb["test_with_lig.pdb"])
    pl_complex.set_ligand("EST:B")
    return pl_complex


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


def test_ligand_and_receptor_indices(pl_complex):
    pl = deepcopy(pl_complex)
    pl._ligand_and_receptor_indices()

    expected_receptor = list(range(0, 146))
    assert len(pl._receptor_indices) == len(expected_receptor)
    assert pl._receptor_indices == expected_receptor

    expected_ligand = list(range(146, 166))
    assert len(pl._lig_indices) == len(expected_ligand)
    assert pl._lig_indices == expected_ligand


def test_ligand_to_mol(pl_complex):
    pl = deepcopy(pl_complex)
    pl._ligand_and_receptor_indices()
    pl._ligand_to_mol()
    assert pl.ligand.GetNumAtoms() == 20


def test_ligand_to_mol_empty_indices_list_raises_error(pl_complex):
    with pytest.raises(exc.NoLigandIndicesError):
        pl_complex._ligand_to_mol()


def test_remove_ligand(pl_complex):
    pl = deepcopy(pl_complex)
    pl._ligand_and_receptor_indices()
    pl._remove_ligand()
    assert pl.traj.n_atoms == 146
    assert pl.traj.n_chains == 1
    assert pl.topology.n_atoms == 146
    assert pl.topology.n_chains == 1


def test_remove_ligand_empty_indices_list_raises_error(pl_complex):
    with pytest.raises(exc.NoLigandIndicesError):
        pl_complex._remove_ligand()


def test_has_hydrogens(pl_complex):
    # TODO: test with a trajectory that contains hydrogen
    assert not pl_complex.has_hydrogens()


def test_add_hydrogens(mocker):
    pass


def test_fix_ligand_smiles_is_given(mocker, pl_complex):
    mock_assign = mocker.patch(
        "openpharmacophore.pharmacophore.pl_complex.Chem.AssignBondOrdersFromTemplate")
    mock_assign.return_value.GetNumAtoms.return_value = 20

    smiles = "CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O"
    pl = deepcopy(pl_complex)
    pl._ligand = mocker.Mock()
    pl._ligand.GetNumAtoms.return_value = 20
    pl.fix_ligand(smiles)

    args = mock_assign.call_args_list
    assert len(args) == 1

    args = args[0][0]
    assert args[0].GetNumAtoms() == 20
    assert pl.ligand.GetNumAtoms() == 20


def test_fix_ligand_no_smiles_given(mocker, pl_complex):
    mock_assign = mocker.patch(
        "openpharmacophore.pharmacophore.pl_complex.Chem.AssignBondOrdersFromTemplate")
    mock_assign.return_value.GetNumAtoms.return_value = 20
    mocker.patch(
        "openpharmacophore.pharmacophore.pl_complex.PLComplex._pdb_id_to_smi",
        return_value="CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O"
    )

    pl = deepcopy(pl_complex)
    pl._ligand = mocker.Mock()
    pl._ligand.GetNumAtoms.return_value = 20
    pl.fix_ligand()

    args = mock_assign.call_args_list
    assert len(args) == 1

    args = args[0][0]
    assert args[0].GetNumAtoms() == 20
    assert pl.ligand.GetNumAtoms() == 20


def test_pdb_id_to_smi(mocker):
    mock_open = mocker.patch(
        "openpharmacophore.pharmacophore.pl_complex.open",
        new=mocker.mock_open())
    mock_open.return_value.readlines.return_value = [
        "EST CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O",
        "DAO CCCCCCCCCCCC(=O)O"
    ]

    smiles = PLComplex._pdb_id_to_smi("EST:B")
    assert smiles == "CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O"
    mock_open.assert_called_once_with(data.pdb_to_smi)


def test_fix_ligand_no_smiles_found(mocker):
    mock_open = mocker.patch(
        "openpharmacophore.pharmacophore.pl_complex.open",
        new=mocker.mock_open())
    mock_open.return_value.readlines.return_value = [
        "EST CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O",
        "DAO CCCCCCCCCCCC(=O)O"
    ]

    pl = PLComplex(data.pdb["test_with_lig.pdb"])
    pl._ligand_ids = ["ATP"]
    pl.set_ligand("ATP")
    with pytest.raises(exc.SmilesNotFoundError):
        pl.fix_ligand()


def test_fix_ligand_pdb_id_unl():
    pl = PLComplex(data.pdb["test_with_lig.pdb"])
    pl._ligand_ids = ["UNL"]
    pl.set_ligand("UNL")

    with pytest.raises(exc.SmilesNotFoundError):
        pl.fix_ligand()


def test_fix_ligand_template_and_lig_atom_number_different():
    pl = PLComplex(data.pdb["test_with_lig.pdb"])
    pl._ligand = Chem.MolFromSmiles("CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O")
    with pytest.raises(exc.DifferentNumAtomsError):
        pl.fix_ligand(smiles="CCCCCCCCCCCC(=O)O")
