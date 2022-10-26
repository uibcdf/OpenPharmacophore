import openpharmacophore as oph
import openpharmacophore.data as data
import nglview as nv
import pytest

pytest.skip("Extract method not implemented yet",
            allow_module_level=True)


def assert_view_contains_pharmacophore(view, n_points):
    assert "nglview.adaptor.MDTrajTrajectory" in view._ngl_component_names
    # There is at least one sphere
    assert "nglview.shape.Shape" in view._ngl_component_names

    n_shapes = len([comp for comp in view._ngl_component_names
                    if comp == "nglview.shape.shape"])
    assert n_shapes >= n_points  # Direction arrows can make n_shapes greater


def test_ligand_receptor_pharmacophore_hydrogen_bonding_points():
    # We create a pharmacophore that contains only hydrogen donor and acceptor points
    pharmacophore = oph.load(data.pdb["3bbh_A_chain.pdb"])
    assert isinstance(pharmacophore, oph.LigandReceptorPharmacophore)
    # We know that the file contains a single ligand
    lig_ids = pharmacophore.receptor.ligand_ids
    assert lig_ids == ["SFG:B"]

    # We extract the pharmacophore
    pharmacophore.extract(lig_ids[0],
                          features=["hb donor", "hb acceptor"])
    assert len(pharmacophore[0]) > 0
    feat_names = [p.feat_name for p in pharmacophore[0]]
    for name in feat_names:
        assert name == "hb acceptor" or name == "hb donor"

    # We inspect the ligand to see that it was correctly extracted.
    # SFG has 50 atoms including hydrogens
    assert pharmacophore.ligand.GetNumAtoms() == 50

    # The protein should contain hydrogens
    assert pharmacophore.receptor.has_hydrogens()

    # Finally we visualize the pharmacophore.
    view = pharmacophore.show(ligand=True, receptor=True)
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))

    # We create a view of the binding site, so we can the residues that
    # are involved in hydrogen bonding clearly.
    view = pharmacophore.bsite_view()
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))


def test_ligand_receptor_pharmacophore_hydrophobic_points():
    # We create a pharmacophore that contains only hydrogen donor and acceptor points
    pharmacophore = oph.load(data.pdb["1m7w_A_chain.pdb"])
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
    assert all([p.feat_name == "hydrophobicity" for p in pharmacophore[0]])

    # We inspect the ligand to see that it was correctly extracted.
    # DAO has 14 atoms excluding hydrogens
    assert pharmacophore.ligand.GetNumAtoms() == 14

    # The protein should not contain hydrogens
    assert not pharmacophore.receptor.has_hydrogens()

    # Finally we visualize the pharmacophore.
    view = pharmacophore.show(ligand=True, receptor=True)
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))

    # We create a view of the binding site, so we can the residues that
    # are involved in hydrophobic interactions clearly.
    view = pharmacophore.bsite_view()
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))


def test_ligand_receptor_pharmacophore_aromatic_points():
    # We create a pharmacophore that contains only aromatic points
    pharmacophore = oph.load(data.pdb["1xdn.pdb"])
    assert isinstance(pharmacophore, oph.LigandReceptorPharmacophore)
    # We know that the file contains a single ligand
    lig_ids = pharmacophore.receptor.ligand_ids
    assert lig_ids == ["ATP:B"]

    # We extract the pharmacophore
    pharmacophore.extract(lig_ids[0],
                          features=["aromatic ring"],
                          add_hydrogens=False)
    assert len(pharmacophore[0]) == 1
    assert all([p.feat_name == "aromatic ring" for p in pharmacophore[0]])

    # We inspect the ligand to see that it was correctly extracted.
    # ATP has 31 atoms excluding hydrogens
    assert pharmacophore.ligand.GetNumAtoms() == 31

    # The protein should not contain hydrogens
    assert not pharmacophore.receptor.has_hydrogens()

    # Finally we add the pharmacophore to an existing view.
    view = nv.NGLWidget()
    pharmacophore.add_to_view(view)
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))

    # We create a view of the binding site, so we can see the residues that
    # are involved in aromatic interactions.
    view = pharmacophore.bsite_view()
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))


def test_ligand_receptor_pharmacophore_from_pdb():
    # We want to create a pharmacophore for the protein-ligand complex of
    # estrogen receptor with estradiol.

    # We load a pharmacophore from a pdb file. This file contains
    # a single peptide chain with a ligand, it contains no hydrogens.
    pharmacophore = oph.load(data.pdb["er_alpha_A_chain.pdb"])
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
    assert pharmacophore.ligand.GetNumAtoms() == 44

    # The protein should contain hydrogens
    assert pharmacophore.receptor.has_hydrogens()

    # Finally we want to view our pharmacophore using nglview.
    view = pharmacophore.show(ligand=True, receptor=True)
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))

    # We create a view of the binding site, so we can see the residues that
    # are involved in protein ligand interactions clearly.
    view = pharmacophore.bsite_view()
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))

# TODO: add a test for a pdb that contains hydrogens and its ligand has correct bonds
