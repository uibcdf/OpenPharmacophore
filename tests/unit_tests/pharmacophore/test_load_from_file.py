from openpharmacophore import load_from_file, LigandReceptorPharmacophore, LigandBasedPharmacophore
import openpharmacophore.data as data


def test_load_ligand_receptor_pharma_ph4():
    pharmacophore = load_from_file(data.pharmacophores["gmp.ph4"])
    assert isinstance(pharmacophore, LigandReceptorPharmacophore)
    assert len(pharmacophore) == 1
    assert len(pharmacophore[0]) == 10


def test_load_ligand_receptor_pharma_mol2():
    pharmacophore = load_from_file(data.pharmacophores["elastase.mol2"])
    assert isinstance(pharmacophore, LigandReceptorPharmacophore)
    assert len(pharmacophore) == 8
    assert pharmacophore.num_frames == 8
    assert pharmacophore._pharmacophores_frames == list(range(0, 8))


def test_load_ligand_receptor_pharma_json():
    pharmacophore = load_from_file(data.pharmacophores["1M70.json"])
    assert isinstance(pharmacophore, LigandReceptorPharmacophore)
    assert len(pharmacophore) == 1
    assert len(pharmacophore[0]) == 5


def test_load_ligand_receptor_pharma_pml():
    pharmacophore = load_from_file(data.pharmacophores["ligscout.pml"])
    assert isinstance(pharmacophore, LigandReceptorPharmacophore)
    assert len(pharmacophore) == 1
    assert len(pharmacophore) == 4


def test_load_ligand_based_pharma_ph4():
    pharmacophore = load_from_file(data.pharmacophores["gmp.ph4"],
                                   pharma_type="ligand")
    assert isinstance(pharmacophore, LigandBasedPharmacophore)
    assert len(pharmacophore) == 1
    assert len(pharmacophore) == 10


def test_load_ligand_based_pharma_mol2():
    pharmacophore = load_from_file(data.pharmacophores["elastase.mol2"],
                                   pharma_type="ligand")
    assert isinstance(pharmacophore, LigandBasedPharmacophore)
    assert len(pharmacophore) == 4


def test_load_ligand_based_pharma_json():
    pharmacophore = load_from_file(data.pharmacophores["1M70.json"],
                                   pharma_type="ligand")
    assert isinstance(pharmacophore, LigandBasedPharmacophore)
    assert len(pharmacophore) == 1
    assert len(pharmacophore[0]) == 5


def test_load_ligand_based_pharma_pml():
    pharmacophore = load_from_file(data.pharmacophores["ligscout.pml"],
                                   pharma_type="ligand")
    assert isinstance(pharmacophore, LigandBasedPharmacophore)
    assert len(pharmacophore) == 1
    assert len(pharmacophore) == 4
