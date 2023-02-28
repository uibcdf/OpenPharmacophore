import openpharmacophore as oph
import pytest


@pytest.mark.skip(reason="Not implemented yet")
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
    assert len(pharmacophore) == 10

    # We inspect the features of a pharmacophore. We expect thrombin pharmacophore
    # to have an aromatic ring and acceptor and a donor
    feat_names = [p.feature_name for p in pharmacophore[0]]
    assert "hb acceptor" in feat_names
    assert "aromatic ring" in feat_names
    assert "hb acceptor" in feat_names

    # We visualize the best scoring pharmacophore
    viewer = oph.Viewer(pharmacophore=pharmacophore, ligands=ligands)
    view = viewer.show()
    # The view should contain a component for the ligand + 3 the pharmacophoric
    # points
    assert len(view._ngl_component_names) == 4
