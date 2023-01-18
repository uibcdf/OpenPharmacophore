import openpharmacophore as oph


def test_binding_site_complex_structure_has_hydrogens(pdb_3bbh_with_hydrogen):
    # We want to extract the binding site of a protein-ligand complex.
    # The complex structure already contains hydrogens, so there is no need to add them.
    protein = oph.load(pdb_3bbh_with_hydrogen)
    assert protein.has_hydrogens
    assert protein.has_ligands

    lig_ids = protein.ligand_ids
    assert len(lig_ids) == 1

    # We set the ligand whose chemical features we want to obtain, and we
    # see that its bond orders are not correct
    ligand = protein.get_ligand(lig_ids[0])
    assert ligand.n_atoms == 50
    assert ligand.has_hydrogens
    assert not ligand.has_aromatic_bonds()

    # We fix the ligand. Its bond orders are correct now
    ligand.fix_bond_order(
        smiles="c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)C[C@H](CC[C@@H](C(=O)O)N)N)O)O)N",
    )
    assert ligand.has_aromatic_bonds()
    assert ligand.get_conformer(0).shape == (50, 3)

    # We create a binding site to obtain the receptor chemical features
    bsite = oph.ComplexBindingSite(protein, ligand)

    # Extract chemical features and visualize them
    receptor_feats = bsite.get_chem_feats(frame=0)
    ligand_feats = ligand.get_chem_feats(conformer=0)
    assert len(receptor_feats) > 0
    assert len(ligand_feats) > 0

    # viewer = oph.Viewer(protein=protein, ligands=ligand)
    # viewer.add_chem_feats([ligand_feats, receptor_feats])
    # viewer.set_protein_indices(bsite_extractor.get_indices(),
    #                            frame=0)
    # viewer.show()


def test_binding_site_complex_structure_does_not_contain_hydrogen(pdb_1m7w):
    # We want to prepare the protein ligand complex of the pdb 1M7W
    # to be ready for pharmacophore extraction. This complex contains
    # lauric acid (id DAO, C12H24O2) as its ligand. This structure does not contain
    # hydrogen, so we need to add them to the ligand and to the receptor.
    protein = oph.load(pdb_1m7w)
    n_atoms_start = protein.n_atoms
    assert protein.has_ligands

    lig_ids = protein.ligand_ids
    assert lig_ids == ["DAO:B"]
    assert not protein.has_hydrogens

    # First we need to extract the ligand, so we can fix its bond orders
    # and add hydrogens to it.
    ligand = protein.get_ligand(lig_ids[0])
    assert ligand.n_atoms == 14
    # We need to remove the ligand from the protein because it does not contain hydrogen
    protein.remove_ligand(lig_ids[0])
    # The protein should have 14 less atoms because the ligand was removed
    assert protein.n_atoms == n_atoms_start - 14
    # The ligand does not have hydrogens yet and all its bonds are single
    assert not ligand.has_hydrogens
    assert not ligand.has_double_bonds()

    # After fixing ligand, it must have a double bond and hydrogens
    ligand.fix_bond_order(smiles="CCCCCCCCCCCC(=O)O")
    assert ligand.has_double_bonds()

    ligand.add_hydrogens()
    assert ligand.n_atoms == 38
    assert ligand.has_hydrogens
    assert ligand.n_conformers == 1
    assert ligand.get_conformer(0).shape == (38, 3)

    # We add hydrogens to the receptor
    protein.add_hydrogens()
    assert protein.has_hydrogens

    # We can create the binding site
    bsite = oph.ComplexBindingSite(protein, ligand)

    # Extract chemical features and visualize them
    receptor_feats = bsite.get_chem_feats(frame=0)
    ligand_feats = ligand.get_chem_feats(conformer=0)
    assert len(receptor_feats) > 0
    assert len(ligand_feats) > 0
