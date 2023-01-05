from openpharmacophore import load, LigandReceptorPharmacophore, LigandBasedPharmacophore


def test_load_ligand_based_pharmacophore_from_mol_file(ligands_smi):
    pharmacophore = load(ligands_smi)
    assert isinstance(pharmacophore, LigandBasedPharmacophore)
    assert len(pharmacophore) == 0
    assert len(pharmacophore.ligands) == 5


def test_load_ligand_receptor_pharmacophore_from_traj_file(pdb_1ncr_path):
    pharmacophore = load(pdb_1ncr_path)
    assert isinstance(pharmacophore, LigandReceptorPharmacophore)
    assert pharmacophore.num_frames == 0
