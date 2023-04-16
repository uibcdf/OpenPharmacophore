from io import StringIO
from rdkit import Chem
from openpharmacophore import MolDB


def test_iterate_smiles_list():
    smiles = [
        "CN=C=O",
        "COc1cc(C=O)ccc1O",
    ]
    mol_db = MolDB()
    mol_db.from_smiles(smiles)

    molecules = []
    for mol in mol_db:
        molecules.append(mol)

    assert len(molecules) == 2
    assert molecules[0].n_atoms == 4
    assert molecules[1].n_atoms == 11


def test_iterate_file():
    file_ = StringIO("CN=C=O\nCOc1cc(C=O)ccc1O")
    mol_db = MolDB()
    mol_db._from_file(file_, extension="smi")

    molecules = []
    for mol in mol_db:
        molecules.append(mol)

    assert len(molecules) == 2
    assert molecules[0].n_atoms == 4
    assert molecules[1].n_atoms == 11


def test_iterate_file_list():
    files = [
        StringIO("CN=C=O\nCOc1cc(C=O)ccc1O"),
        StringIO("CN1CCC[C@H]1c2cccnc2")
    ]

    mol_db = MolDB()
    mol_db._from_file_list(files)

    molecules = []
    for mol in mol_db:
        molecules.append(mol)

    assert len(molecules) == 3
    assert molecules[0].n_atoms == 4
    assert molecules[1].n_atoms == 11
    assert molecules[2].n_atoms == 12


def test_iterate_rdkit_mol():
    rdkit_mols = [
        Chem.MolFromSmiles("CN=C=O"),
        Chem.MolFromSmiles("COc1cc(C=O)ccc1O"),
    ]

    mol_db = MolDB()
    mol_db.from_rdkit(rdkit_mols)

    molecules = []
    for mol in mol_db:
        molecules.append(mol)

    assert len(molecules) == 2
    assert molecules[0].n_atoms == 4
    assert molecules[1].n_atoms == 11
