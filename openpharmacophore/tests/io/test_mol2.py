import openpharmacophore.data as data
import openpharmacophore.io as io


def test_load_mol2_file():
    file_name = data.ligands["ace"]
    molecules = io.load_mol2_file(file_name=file_name)

    assert len(molecules) == 3
    assert molecules[0].GetNumAtoms() == 14
    assert molecules[1].GetNumAtoms() == 25
    assert molecules[2].GetNumAtoms() == 29
