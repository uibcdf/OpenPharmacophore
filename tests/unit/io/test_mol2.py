from openpharmacophore.io.mol2 import load_mol2_ligands, mol2_iterator


def test_load_mol2_file(ligands_mol2):
    molecules = load_mol2_ligands(file_name=ligands_mol2)

    assert len(molecules) == 3
    assert molecules[0].GetNumAtoms() == 14
    assert molecules[1].GetNumAtoms() == 25
    assert molecules[2].GetNumAtoms() == 29


def test_mol2_iterator(ligands_mol2):
    molecules = []
    for mol in mol2_iterator(ligands_mol2):
        molecules.append(mol)

    assert len(molecules) == 3
    assert molecules[0].GetNumAtoms() == 14
    assert molecules[1].GetNumAtoms() == 25
    assert molecules[2].GetNumAtoms() == 29
