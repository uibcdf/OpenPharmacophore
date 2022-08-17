import openpharmacophore.io as io
import openpharmacophore.data as data
import pytest
from rdkit import Chem


@pytest.mark.parametrize("file_name", [
    "ace.mol2",
    "clique_detection.smi"
])
def test_load_molecules_file(file_name):
    ligands_ = io.load_molecules_file(file_name=data.ligands["clique_detection"],
                                      titleLine=False)
    assert len(ligands_) == 5

    ligands_ = io.load_molecules_file(file_name=data.ligands["ace"])
    assert len(ligands_) == 3

    assert isinstance(ligands_, list)
    for lig in ligands_:
        assert isinstance(lig, Chem.Mol)


def test_smi_has_header_and_id():
    file = data.ligands["mols"]
    has_header, has_id = io.smi_has_header_and_id(file)
    assert not has_header
    assert not has_id

    file = data.ligands["BAAAML"]
    has_header, has_id = io.smi_has_header_and_id(file)
    assert has_header
    assert has_id


def test_mol2_mol_generator():
    mol2_file = data.ligands["ace"]

    molecules = []
    with open(mol2_file, "r") as fp:
        for mol in io.mol2_mol_generator(fp):
            molecules.append(mol)

    assert len(molecules) == 3
    assert molecules[0].GetNumAtoms() == 14
    assert molecules[1].GetNumAtoms() == 25
    assert molecules[2].GetNumAtoms() == 29


def test_smiles_mol_generator():
    smi_file = data.ligands["mols"]

    n_molecules = 0
    with open(smi_file, "r") as fp:
        for mol in io.smiles_mol_generator(fp, header=False, mol_id=False):
            n_molecules += 1
            assert isinstance(mol, Chem.Mol)

    assert n_molecules == 5

    smi_file = data.ligands["BAAAML"]

    n_molecules = 0
    with open(smi_file, "r") as fp:
        for mol in io.smiles_mol_generator(fp, header=True, mol_id=True):
            n_molecules += 1
            assert isinstance(mol, Chem.Mol)

    assert n_molecules == 23
