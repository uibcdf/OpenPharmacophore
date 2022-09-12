import openpharmacophore.data as data


def test_all_data_dictionaries_are_populated():

    assert len(data.ligands) == 3
    assert len(data.pdb) == 5
    assert len(data.pharmacophores) == 5
    # assert len(data.trajectories) == 1
    # assert len(data.topologies) == 1


def test_can_load_a_file_from_each_dictionary():

    molecules = data.ligands["mols"]
    with open(molecules) as fp:
        lines = fp.readlines()
        assert len(lines) == 5

    pdb_1ncr = data.pdb["1ncr"]
    with open(pdb_1ncr) as fp:
        header = fp.readline()
        assert "HEADER" in header

    pharmacophore = data.pharmacophores["elastase"]
    with open(pharmacophore) as fp:
        header = fp.readline()
        assert "@<TRIPOS>" in header
