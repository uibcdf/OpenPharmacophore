import numpy as np
import openpharmacophore._private_tools.exceptions as exc
import openpharmacophore.data as data
from openpharmacophore.pharmacophore.pl_complex import PLComplex
import pytest
from copy import deepcopy
from rdkit import Chem
from openmm.app import PDBFile


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
    pl.ligand_and_receptor_indices()

    expected_receptor = list(range(0, 146))
    assert len(pl._receptor_indices) == len(expected_receptor)
    assert pl._receptor_indices == expected_receptor

    expected_ligand = list(range(146, 166))
    assert len(pl._lig_indices) == len(expected_ligand)
    assert pl._lig_indices == expected_ligand


def test_ligand_to_mol(pl_complex):
    pl = deepcopy(pl_complex)
    pl.ligand_and_receptor_indices()
    pl.ligand_to_mol()
    assert pl.ligand.GetNumAtoms() == 20


def test_ligand_to_mol_empty_indices_list_raises_error(pl_complex):
    with pytest.raises(exc.NoLigandIndicesError):
        pl_complex.ligand_to_mol()


def test_remove_ligand(pl_complex):
    pl = deepcopy(pl_complex)
    pl.ligand_and_receptor_indices()
    pl.remove_ligand()
    assert pl.traj.n_atoms == 146
    assert pl.traj.n_chains == 1
    assert pl.topology.n_atoms == 146
    assert pl.topology.n_chains == 1


def test_remove_ligand_empty_indices_list_raises_error(pl_complex):
    with pytest.raises(exc.NoLigandIndicesError):
        pl_complex.remove_ligand()


def test_has_hydrogens(pl_complex):
    # TODO: test with a trajectory that contains hydrogen
    assert not pl_complex.has_hydrogens()


def test_add_hydrogens(mocker, pl_complex):
    mock_model_to_traj = mocker.patch(
        "openpharmacophore.pharmacophore.pl_complex.PLComplex._modeller_to_trajectory"
    )
    mock_modeller = mocker.patch(
        "openpharmacophore.pharmacophore.pl_complex.Modeller"
    )
    pl_complex.add_hydrogens()

    mock_modeller.assert_called_once()
    _, _, kwargs = mock_modeller.mock_calls[0]
    assert len(kwargs["positions"]) == 166
    assert kwargs["topology"].getNumAtoms() == 166

    mock_model_to_traj.assert_called_once_with(mock_modeller.return_value)


@pytest.fixture()
def estradiol_mol():
    """ Returns estradiol as an rdkit.Mol object.
    """
    pdb_string = """
REMARK   1 CREATED WITH MDTraj 1.9.7, 2022-10-17
CRYST1  105.500  105.500  136.080  90.00  90.00 120.00 P 1           1 
MODEL        0
ATOM   5944  C1  EST A 600     104.106  17.203  24.775  1.00  0.00           C  
ATOM   5945  C2  EST A 600     102.995  17.834  25.370  1.00  0.00           C  
ATOM   5946  C3  EST A 600     101.695  17.355  25.120  1.00  0.00           C  
ATOM   5947  O3  EST A 600     100.598  17.990  25.704  1.00  0.00           O  
ATOM   5948  C4  EST A 600     101.506  16.240  24.274  1.00  0.00           C  
ATOM   5949  C5  EST A 600     102.621  15.588  23.660  1.00  0.00           C  
ATOM   5950  C6  EST A 600     102.371  14.379  22.735  1.00  0.00           C  
ATOM   5951  C7  EST A 600     103.644  13.753  22.086  1.00  0.00           C  
ATOM   5952  C8  EST A 600     104.898  13.873  22.953  1.00  0.00           C  
ATOM   5953  C9  EST A 600     105.178  15.388  23.261  1.00  0.00           C  
ATOM   5954  C10 EST A 600     103.957  16.078  23.918  1.00  0.00           C  
ATOM   5955  C11 EST A 600     106.462  15.459  24.125  1.00  0.00           C  
ATOM   5956  C12 EST A 600     107.711  14.803  23.508  1.00  0.00           C  
ATOM   5957  C13 EST A 600     107.463  13.343  23.124  1.00  0.00           C  
ATOM   5958  C14 EST A 600     106.170  13.270  22.242  1.00  0.00           C  
ATOM   5959  C15 EST A 600     106.228  11.821  21.792  1.00  0.00           C  
ATOM   5960  C16 EST A 600     107.701  11.713  21.263  1.00  0.00           C  
ATOM   5961  C17 EST A 600     108.494  12.719  22.135  1.00  0.00           C  
ATOM   5962  O17 EST A 600     109.610  12.027  22.746  1.00  0.00           O  
ATOM   5963  C18 EST A 600     107.379  12.449  24.419  1.00  0.00           C  
TER    5964      EST A 600
ENDMDL
CONECT    1    2   11
CONECT    2    1    3
CONECT    3    2    4    5
CONECT    4    3
CONECT    5    3    6
CONECT    6    5    7   11
CONECT    7    6    8
CONECT    8    7    9
CONECT    9    8   10   15
CONECT   10    9   11   12
CONECT   11    1    6   10
CONECT   12   10   13
CONECT   13   12   14
CONECT   14   13   15   18   20
CONECT   15    9   14   16
CONECT   16   15   17
CONECT   17   16   18
CONECT   18   14   17   19
CONECT   19   18
CONECT   20   14
END
    """
    estradiol = Chem.MolFromPDBBlock(pdb_string)
    assert estradiol is not None
    return estradiol


def test_fix_ligand_smiles_is_given(pl_complex, estradiol_mol):

    smiles = "CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O"
    pl = deepcopy(pl_complex)
    pl._ligand = deepcopy(estradiol_mol)
    assert not all(b.GetIsAromatic() for b in pl.ligand.GetBonds())
    assert len([a.GetSymbol() for a in pl.ligand.GetAtoms()
                if a.GetSymbol() == "H"]) == 0

    pl.fix_ligand(smiles=smiles)
    assert pl.ligand.GetNumAtoms() == 44
    assert any(b.GetIsAromatic() for b in pl.ligand.GetBonds())
    assert len([a.GetSymbol() for a in pl.ligand.GetAtoms()
                if a.GetSymbol() == "H"]) == 24


def test_fix_ligand_no_smiles_given(mocker, pl_complex, estradiol_mol):
    mocker.patch(
        "openpharmacophore.pharmacophore.pl_complex.PLComplex._pdb_id_to_smi",
        return_value="C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O"
    )

    pl = deepcopy(pl_complex)
    pl._ligand = deepcopy(estradiol_mol)
    assert not all(b.GetIsAromatic() for b in pl.ligand.GetBonds())
    assert len([a.GetSymbol() for a in pl.ligand.GetAtoms()
                if a.GetSymbol() == "H"]) == 0

    pl.fix_ligand()
    assert pl.ligand.GetNumAtoms() == 44
    assert any(b.GetIsAromatic() for b in pl.ligand.GetBonds())
    assert len([a.GetSymbol() for a in pl.ligand.GetAtoms()
                if a.GetSymbol() == "H"]) == 24


def test_pdb_id_to_smi(mocker):
    mock_open = mocker.patch(
        "openpharmacophore.pharmacophore.pl_complex.open",
        new=mocker.mock_open())
    mock_open.return_value.readlines.return_value = [
        "EST C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O",
        "DAO CCCCCCCCCCCC(=O)O"
    ]

    smiles = PLComplex._pdb_id_to_smi("EST:B")
    assert smiles == "C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O"
    mock_open.assert_called_once_with(data.pdb_to_smi)


def test_fix_ligand_no_smiles_found(mocker):
    mock_open = mocker.patch(
        "openpharmacophore.pharmacophore.pl_complex.open",
        new=mocker.mock_open())
    mock_open.return_value.readlines.return_value = [
        "EST C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O",
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
    pl._ligand = Chem.MolFromSmiles("C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O")
    with pytest.raises(exc.DifferentNumAtomsError):
        pl.fix_ligand(smiles="CCCCCCCCCCCC(=O)O")


def test_modeller_to_trajectory():
    modeller = PDBFile(data.pdb["test_no_lig.pdb"])
    traj = PLComplex._modeller_to_trajectory(modeller)
    assert traj.n_frames == 1
    assert traj.n_atoms == 19
    assert traj.n_chains == 1
    assert traj.topology.n_atoms == 19
    assert traj.topology.n_chains == 1
    assert traj.topology.n_residues == 2

    expected_coords = np.array([
        [44.235,  80.308,  18.419],
        [43.549,  79.243,  17.706],
        [44.528,  78.252,  17.077],
        [45.699,  78.559,  16.853],
        [42.608,  79.792,  16.611],
        [43.375,  80.468,  15.608],
        [41.586,  80.762,  17.208],
        [44.030,  77.052,  16.814],
        [44.799,  76.021,  16.156],
        [44.189,  75.782,  14.791],
        [43.007,  76.033,  14.583],
        [44.721,  74.725,  16.954],
        [45.253,  74.836,  18.355],
        [46.598,  74.629,  18.621],
        [44.412,  75.147,  19.416],
        [47.098,  74.729,  19.906],
        [44.895,  75.245,  20.707],
        [46.246,  75.036,  20.946],
        [46.748,  75.126,  22.224],
    ]) / 10  # convert to nanometers

    assert np.allclose(traj.xyz, expected_coords)


def test_add_fixed_ligand():
    assert False, "Complete me!"
