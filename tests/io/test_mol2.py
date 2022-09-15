import openpharmacophore.data as data
from openpharmacophore.io.mol2 import load_mol2_ligands


def test_load_mol2_file():
    file_name = data.ligands["ace"]
    molecules = load_mol2_ligands(file_name=file_name)

    assert len(molecules) == 3
    assert molecules[0].GetNumAtoms() == 14
    assert molecules[1].GetNumAtoms() == 25
    assert molecules[2].GetNumAtoms() == 29
