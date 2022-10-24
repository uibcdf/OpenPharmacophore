import openpharmacophore.data as data


def test_all_data_dictionaries_are_populated():

    assert "ace.mol2" in data.ligands
    assert "1ncr.pdb" in data.pdb
    assert "elastase.mol2" in data.pharmacophores
    assert "tranches_2D.txt" in data.zinc
    assert "pdb_to_smi.txt" in data.pdb_to_smi
