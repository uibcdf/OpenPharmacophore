from openpharmacophore.io import mol_file_to_list
import openpharmacophore.data as data
from rdkit import Chem


def assert_is_mol_list(mol_list):
    for mol in mol_list:
        assert isinstance(mol, Chem.Mol)


def test_mol_file_to_list_smi():
    molecules = mol_file_to_list(data.ligands["clique_detection"])
    assert len(molecules) == 5
    assert_is_mol_list(molecules)


def test_mol_file_to_list_mol2():
    molecules = mol_file_to_list(data.ligands["ace"])
    assert len(molecules) == 3
    assert_is_mol_list(molecules)


def test_mol_file_to_list_sdf():
    molecules = mol_file_to_list(data.ligands["sdf_example"])
    assert len(molecules) == 3
    assert_is_mol_list(molecules)
