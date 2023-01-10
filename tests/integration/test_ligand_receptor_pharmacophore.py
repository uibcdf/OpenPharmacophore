import openpharmacophore as oph
from assert_view import assert_view_contains_pharmacophore


def test_ligand_receptor_pharmacophore_from_pdb(pdb_er_alpha):
    # We want to create a pharmacophore for the protein-ligand complex of
    # estrogen receptor with estradiol.
    protein = oph.load(pdb_er_alpha)
    # We know that the file contains a single ligand
    assert protein.has_ligands
    lig_ids = protein.ligand_ids
    assert lig_ids == ["EST:B"]

    # We obtain the smiles of the ligand. Necessary to fix its bond order later
    smiles = oph.smiles_from_pdb_id(lig_ids[0])

    # We extract the ligand and fix its bond order and add hydrogens
    ligand = protein.get_ligand(lig_ids[0])
    ligand.fix_bond_order(smiles=smiles)
    ligand.add_hydrogens()
    # Estradiol has 18 Carbons, 2 Oxygens and 24 Hydrogens.
    assert ligand.has_aromatic_bonds()
    assert ligand.n_atoms == 44  # This ligand should contain 50 atoms

    # We add hydrogens to the protein
    protein.add_hydrogens()
    assert protein.has_hydrogens

    # We get the binding site
    binding_site = protein.extract_binding_site(ligand=ligand)

    # We extract the pharmacophore
    pharmacophore = oph.LigandReceptorPharmacophore(binding_site, ligand)
    pharmacophore.extract()
    assert len(pharmacophore[0]) > 0
    # We know that the pharmacophore of ERalpha should contain one aromatic point and
    # at least one hydrophobic point and a donor
    feats = [p.feature_name for p in pharmacophore[0]]
    assert "aromatic ring" in feats
    assert "hydrophobicity" in feats
    assert "hb donor" in feats

    # Finally we visualize the pharmacophore.
    viewer = oph.Viewer(pharmacophore=pharmacophore,
                        receptor=protein,
                        ligands=ligand
                        )
    view = viewer.show(pharma_index=0)
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))

    # We create a view of the binding site, so we can the residues that
    # are involved in protein ligand-interactions.
    view.remove_protein()
    view.add_bsite(binding_site)
    assert_view_contains_pharmacophore(view, len(pharmacophore[0]))
