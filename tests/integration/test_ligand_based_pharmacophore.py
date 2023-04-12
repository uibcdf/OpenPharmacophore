import openpharmacophore as oph


def test_ligand_based_pharmacophore_extraction(thrombin_ligands):
    # We want to extract a ligand based pharmacophore for thrombin.
    # We start by loading the ligands from a sdf file. The ligands in this file
    # are already prepared for extraction with multiple conformers and hydrogens.
    ligands = oph.load(thrombin_ligands)
    assert len(ligands) == 7

    pharmacophore = oph.LigandBasedPharmacophore(ligands)
    # We extract pharmacophores of 3 points and visualize them
    pharmacophore.extract(
        n_points=3, min_actives=len(ligands), max_pharmacophores=10
    )
    assert len(pharmacophore) > 1

    # We inspect the features of a pharmacophore. We expect thrombin pharmacophore
    # to have an aromatic ring and acceptor and a donor
    feat_names = [p.feature_name for p in pharmacophore[0]]
    assert "hb acceptor" in feat_names
    assert "aromatic ring" in feat_names
    assert "hb acceptor" in feat_names

    # We visualize the best scoring pharmacophore
    top = pharmacophore[0]
    viewer = oph.Viewer()
    viewer.add_components([
        top,
        ligands[top.ref_mol]])
    view = viewer.show(struct=top.ref_mol)
    # The visualization should contain a component for the ligand + 3 the pharmacophoric
    # points
    assert len(view._ngl_component_names) == 4
