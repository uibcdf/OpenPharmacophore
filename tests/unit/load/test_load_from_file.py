from openpharmacophore import load_from_file, LigandReceptorPharmacophore, LigandBasedPharmacophore, Pharmacophore


def test_load_ligand_receptor_pharma_ph4(moe_pharmacophore_path):
    pharmacophore = load_from_file(moe_pharmacophore_path)
    assert isinstance(pharmacophore, LigandReceptorPharmacophore)
    assert len(pharmacophore) == 1
    assert isinstance(pharmacophore[0], Pharmacophore)
    assert len(pharmacophore[0]) == 10


def test_load_ligand_receptor_pharma_mol2(mol2_pharmacophore_path_elastase):
    pharmacophore = load_from_file(mol2_pharmacophore_path_elastase)
    assert isinstance(pharmacophore, LigandReceptorPharmacophore)
    assert len(pharmacophore) == 8
    assert isinstance(pharmacophore[0], Pharmacophore)
    assert pharmacophore.num_frames == 8

    expected_frames = list(range(8))
    frames = [p.ref_struct for p in pharmacophore]
    assert expected_frames == frames


def test_load_ligand_receptor_pharma_json(json_pharmacophore_path):
    pharmacophore = load_from_file(json_pharmacophore_path)
    assert isinstance(pharmacophore, LigandReceptorPharmacophore)
    assert len(pharmacophore) == 1
    assert isinstance(pharmacophore[0], Pharmacophore)
    assert len(pharmacophore[0]) == 5


def test_load_ligand_receptor_pharma_pml(ligand_scout_pharmacophore_path):
    pharmacophore = load_from_file(ligand_scout_pharmacophore_path)
    assert isinstance(pharmacophore, LigandReceptorPharmacophore)
    assert len(pharmacophore) == 1
    assert isinstance(pharmacophore[0], Pharmacophore)
    assert len(pharmacophore[0]) == 4


def test_load_ligand_based_pharma_ph4(moe_pharmacophore_path):
    pharmacophore = load_from_file(moe_pharmacophore_path,
                                   pharma_type="ligand")
    assert isinstance(pharmacophore, LigandBasedPharmacophore)
    assert isinstance(pharmacophore[0], Pharmacophore)
    assert len(pharmacophore[0]) == 10


def test_load_ligand_based_pharma_mol2(mol2_pharmacophore_path_elastase):
    pharmacophore = load_from_file(mol2_pharmacophore_path_elastase,
                                   pharma_type="ligand")
    assert isinstance(pharmacophore, LigandBasedPharmacophore)
    assert isinstance(pharmacophore[0], Pharmacophore)
    assert len(pharmacophore[0]) == 4


def test_load_ligand_based_pharma_json(json_pharmacophore_path):
    pharmacophore = load_from_file(json_pharmacophore_path,
                                   pharma_type="ligand")
    assert isinstance(pharmacophore, LigandBasedPharmacophore)
    assert isinstance(pharmacophore[0], Pharmacophore)
    assert len(pharmacophore[0]) == 5


def test_load_ligand_based_pharma_pml(ligand_scout_pharmacophore_path):
    pharmacophore = load_from_file(ligand_scout_pharmacophore_path,
                                   pharma_type="ligand")
    assert isinstance(pharmacophore, LigandBasedPharmacophore)
    assert isinstance(pharmacophore[0], Pharmacophore)
    assert len(pharmacophore[0]) == 4
