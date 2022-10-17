from openpharmacophore import load, LigandReceptorPharmacophore, LigandBasedPharmacophore
import openpharmacophore.data as data


def test_load_ligand_based_pharmacophore():
    pharmacophore = load(data.ligands["mols.smi"])
    assert isinstance(pharmacophore, LigandBasedPharmacophore)
    assert len(pharmacophore) == 0
    assert len(pharmacophore.ligands) == 5


def test_load_ligand_receptor_pharmacophore():
    pharmacophore = load(data.pdb["1ncr.pdb"])
    assert isinstance(pharmacophore, LigandReceptorPharmacophore)
    assert pharmacophore.num_frames == 1
    