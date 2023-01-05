from openpharmacophore.io import mol_file_to_list, mol_file_iterator
from rdkit import Chem


def assert_is_mol_list(mol_list):
    for mol in mol_list:
        assert isinstance(mol, Chem.Mol)


def test_mol_file_to_list_smi(ligands_smi):
    molecules = mol_file_to_list(ligands_smi)
    assert len(molecules) == 5
    assert_is_mol_list(molecules)


def test_mol_file_to_list_mol2(ligands_mol2):
    molecules = mol_file_to_list(ligands_mol2)
    assert len(molecules) == 3
    assert_is_mol_list(molecules)


def test_mol_file_to_list_sdf(ligands_sdf):
    molecules = mol_file_to_list(ligands_sdf)
    assert len(molecules) == 3
    assert_is_mol_list(molecules)


def assert_iterable_contains_mol(mol_list, n_mols):
    count = 0
    for mol in mol_list:
        count += 1
        assert isinstance(mol, Chem.Mol)
    assert count == n_mols


def test_mol_file_iterator_smi(ligands_smi):
    mol_iterator = mol_file_iterator(ligands_smi)
    assert_iterable_contains_mol(mol_iterator, 5)


def test_mol_file_iterator_mol2(ligands_mol2):
    mol_iterator = mol_file_iterator(ligands_mol2)
    assert_iterable_contains_mol(mol_iterator, 3)


def test_mol_file_iterator_sdf(ligands_sdf):
    mol_iterator = mol_file_iterator(ligands_sdf)
    assert_iterable_contains_mol(mol_iterator, 3)
