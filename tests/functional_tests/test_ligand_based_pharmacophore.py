import openpharmacophore as oph


def test_ligand_based_pharmacophore_from_ligand_set():
    # We want to create a ligand-based pharmacophore for the ligands of Thrombin

    # We start by creating a pharmacophore from a list of smiles representing
    # Thrombin ligands
    ligands_smi = [
        "[H]/N=C(\C1CCC(CC1)CNC(=O)[C@@H]2C=C(CN3N2C(=O)N(C3=O)CC(c4ccccc4)c5ccccc5)C)/N",
        "CN[C@H](Cc1ccccc1)C(=O)N2CCC[C@H]2C(=O)NCC3CCC(CC3)N",
        "c1ccc(cc1)S(=O)(=O)CCN2C(=O)N3CC=C[C@H](N3C2=O)C(=O)NC4CCC(CC4)c5cnc([nH]5)N",
        "c1ccc(cc1)S(=O)(=O)CCN2C(=O)N3CC=C[C@H](N3C2=O)C(=O)NCC4CCC(CC4)N",
        "[H]/N=C(/c1ccc(cc1)C[C@H](C(=O)N2CCCCC2)NC(=O)CNS(=O)(=O)c3ccc4ccccc4c3)\\N",
        "[H]/N=C(\c1ccc2c(c1)cc([nH]2)C(=O)N3CCC(CC3)Cc4ccccc4)/N",
        "CCC1CCN(C(=O)[C@H](CCCNC(N)=[NH2+])NS(=O)(=O)c2cccc3c(N(C)C)cccc23)CC1",
    ]
    pharmacophore = oph.LigandBasedPharmacophore()
    pharmacophore.load_ligands_from_smi(ligands_smi)

    # We confirm that the molecules were correctly loaded
    assert len(pharmacophore.ligands) == 7
    assert all([mol is not None for mol in pharmacophore.ligands])

    # We proceed to add hydrogens to the ligands and generate conformers
    pharmacophore.add_hydrogens(ligands="all")
    pharmacophore.generate_conformers(ligands="all", n_confs=5)

    # We extract pharmacophores of 3 points and visualize it
    pharmacophore.extract(
        n_points=4, min_actives=7
    )
    assert len(pharmacophore) > 0
    assert len(pharmacophore.scores) > 0

    # We inspect the features of the pharmacophore. We expect thrombin pharmacophore
    # to have an aromatic ring, at leas one hydrophobic and a positive charge
    feat_names = [p.feature_name for p in pharmacophore[0]]
    assert "hydrophobicity" in feat_names
    assert "aromatic ring" in feat_names
    assert "positive charge" in feat_names

    # We visualize the best scoring pharmacophore
    view = pharmacophore.show(index=0, ligands=True)
    # The view should contain a component for each ligand (5) + the pharmacophoric
    # points
    assert len(view._ngl_component_names) > 5
