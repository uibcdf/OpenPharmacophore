import openpharmacophore as oph
import pytest


@pytest.mark.skip(reason="Ligand class not implemented yet")
def test_ligand_preparation():
    # We want to prepare the ligands of Thrombin for pharmacophore extraction

    # We start by creating a LigandSetObject object from a list of smiles representing
    # Thrombin ligands
    ligands_smi = [
        r"[H]/N=C(\C1CCC(CC1)CNC(=O)[C@@H]2C=C(CN3N2C(=O)N(C3=O)CC(c4ccccc4)c5ccccc5)C)/N",
        r"CN[C@H](Cc1ccccc1)C(=O)N2CCC[C@H]2C(=O)NCC3CCC(CC3)N",
        r"c1ccc(cc1)S(=O)(=O)CCN2C(=O)N3CC=C[C@H](N3C2=O)C(=O)NC4CCC(CC4)c5cnc([nH]5)N",
        r"c1ccc(cc1)S(=O)(=O)CCN2C(=O)N3CC=C[C@H](N3C2=O)C(=O)NCC4CCC(CC4)N",
    ]
    ligands = oph.load_ligands(ligands_smi, form="smi")
    assert len(ligands) == 4

    # We proceed to add hydrogens to the ligands
    ligands.add_hydrogens(indices="all")
    assert all([lig.has_hydrogen() for lig in ligands])

    # We generate conformers for each
    ligands.generate_conformers(indices="all", n_confs=1)
    assert all([lig.n_confs == 1 for lig in ligands])

    # We draw the ligands
    ligands.draw(n_per_row=4)

    # Now we want to find chemical features in the ligands
    chem_feats = oph.extract_chem_feats(ligands)
    # All ligands have aromatic rings, hydrophobic areas, acceptors and donors
    assert len(chem_feats) == 4
    assert all(["aromatic ring" in f for f in chem_feats])
    assert all(["hb acceptor" in f for f in chem_feats])
    assert all(["hb donor" in f for f in chem_feats])
    assert all(["hydrophobicity" in f for f in chem_feats])
    # Finally we create a 2D representation of the ligands with their
    # chemical features highlighted
    ligands.draw(sub_image_size=(300, 280), chem_feats=chem_feats)
