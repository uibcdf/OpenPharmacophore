import openpharmacophore as oph
from assert_view import assert_view_contains_pharmacophore


def test_ligand_receptor_pharmacophore_hydrogen_bonding_points(pdb_3bbh):
    # We create a pharmacophore that contains only hydrogen donor and acceptor points
    pharmacophore = oph.load(pdb_3bbh)
    assert isinstance(pharmacophore, oph.LigandReceptorPharmacophore)
    # We know that the file contains a single ligand
    lig_ids = pharmacophore.receptor.ligand_ids
    assert lig_ids == ["SFG:B"]

    # We extract the pharmacophore
    pharmacophore.extract(lig_ids[0],
                          features=["hb donor", "hb acceptor"])
    assert len(pharmacophore[0]) > 0
    feat_names = [p.feature_name for p in pharmacophore[0]]
    for name in feat_names:
        assert name == "hb acceptor" or name == "hb donor"

    # We inspect the ligand to see that it was correctly extracted.
    # SFG has 50 atoms including hydrogens
    assert pharmacophore.receptor.ligand.GetNumAtoms() == 50

    # The protein should contain hydrogens
    assert pharmacophore.receptor.has_hydrogens()

    # Finally we visualize the pharmacophore.
    view = pharmacophore.show(ligand=True, receptor=True)
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))

    # We create a view of the binding site, so we can the residues that
    # are involved in hydrogen bonding clearly.
    bsite_indices = pharmacophore.receptor.binding_site_indices(0)
    view = pharmacophore.show(indices=bsite_indices)
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))


def test_ligand_receptor_pharmacophore_hydrophobic_points(pdb_1m7w):
    # We create a pharmacophore that contains only hydrophobic points
    pharmacophore = oph.load(pdb_1m7w)
    assert isinstance(pharmacophore, oph.LigandReceptorPharmacophore)
    # We know that the file contains a single ligand
    lig_ids = pharmacophore.receptor.ligand_ids
    assert lig_ids == ["DAO:B"]

    # We know the smiles of the ligand do we pass it to the extract method
    # so, it can fix the ligand bonds and obtain an accurate pharmacophore
    smiles = "CCCCCCCCCCCC(=O)O"
    pharmacophore.extract(lig_ids[0],
                          features=["hydrophobicity"],
                          smiles=smiles,
                          add_hydrogens=False)
    assert len(pharmacophore[0]) > 0
    assert all([p.feature_name == "hydrophobicity" for p in pharmacophore[0]])

    # We inspect the ligand to see that it was correctly extracted.
    # DAO has 14 atoms excluding hydrogens
    assert pharmacophore.receptor.ligand.GetNumAtoms() == 14

    # The protein should not contain hydrogens
    assert not pharmacophore.receptor.has_hydrogens()

    # Finally we visualize the pharmacophore.
    view = pharmacophore.show(ligand=True, receptor=True)
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))

    # We create a view of the binding site, so we can the residues that
    # are involved in hydrophobic interactions clearly.
    bsite_indices = pharmacophore.receptor.binding_site_indices(0)
    view = pharmacophore.show(indices=bsite_indices)
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))


def test_ligand_receptor_pharmacophore_aromatic_points(pdb_1xdn):
    # We create a pharmacophore that contains only aromatic points
    pharmacophore = oph.load(pdb_1xdn)
    assert isinstance(pharmacophore, oph.LigandReceptorPharmacophore)
    # We know that the file contains a single ligand
    lig_ids = pharmacophore.receptor.ligand_ids
    assert lig_ids == ["ATP:B"]

    # We extract the pharmacophore
    pharmacophore.extract(lig_ids[0],
                          features=["aromatic ring"],
                          add_hydrogens=False)
    assert len(pharmacophore[0]) == 1
    assert all([p.feature_name == "aromatic ring" for p in pharmacophore[0]])

    # We inspect the ligand to see that it was correctly extracted.
    # ATP has 31 atoms excluding hydrogens
    assert pharmacophore.receptor.ligand.GetNumAtoms() == 31

    # The protein should not contain hydrogens
    assert not pharmacophore.receptor.has_hydrogens()

    # Finally we visualize the pharmacophore.
    view = pharmacophore.show()
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))

    # We create a view of the binding site, so we can see the residues that
    # are involved in aromatic interactions.
    bsite_indices = pharmacophore.receptor.binding_site_indices(0)
    view = pharmacophore.show(indices=bsite_indices)
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))


def test_ligand_receptor_pharmacophore_from_pdb(pdb_er_alpha):
    # We want to create a pharmacophore for the protein-ligand complex of
    # estrogen receptor with estradiol.

    # We load a pharmacophore from a pdb file. This file contains
    # a single peptide chain with a ligand, it contains no hydrogens.
    pharmacophore = oph.load(pdb_er_alpha)
    assert isinstance(pharmacophore, oph.LigandReceptorPharmacophore)
    # We call find ligands method to ensure our pdb contains a single ligand.
    # We know estradiol has the id EST.
    lig_ids = pharmacophore.receptor.ligand_ids
    assert lig_ids == ["EST:B"]

    # With the full ligand id we can now extract a ligand-receptor based pharmacophore
    ligand_id = lig_ids[0]
    pharmacophore.extract(ligand_id)
    assert len(pharmacophore[0]) > 0

    # We inspect the ligand to see that it was correctly extracted.
    # Estradiol has 18 Carbons, 2 Oxygens and 24 Hydrogens.
    assert pharmacophore.receptor.ligand.GetNumAtoms() == 44

    # The protein should contain hydrogens
    assert pharmacophore.receptor.has_hydrogens()

    # We know that the pharmacophore of ERalpha should contain one aromatic point and
    # at least one hydrophobic point
    feats = [p.feature_name for p in pharmacophore[0]]
    assert "aromatic ring" in feats
    assert "hydrophobicity" in feats

    # Finally we want to view our pharmacophore using nglview.
    view = pharmacophore.show(ligand=True, receptor=True)
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))

    # We create a view of the binding site, so we can see the residues that
    # are involved in protein ligand interactions clearly.
    bsite_indices = pharmacophore.receptor.binding_site_indices(0)
    view = pharmacophore.show(indices=bsite_indices)
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))
