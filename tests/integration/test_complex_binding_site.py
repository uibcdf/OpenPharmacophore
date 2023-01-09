import openpharmacophore as oph


def test_pl_complex_bsite(pdb_3bbh_with_hydrogen):
    # We want to extract the binding site of a protein-ligand complex
    protein = oph.load(pdb_3bbh_with_hydrogen)
    assert protein.has_hydrogens
    assert protein.has_ligands

    lig_ids = protein.ligand_ids()
    assert len(lig_ids) == 1

    # We set the ligand whose chemical features we want to obtain, and we
    # see that its bond orders are not correct
    ligand = protein.extract(lig_ids[0])
    assert ligand.n_atoms == 27

    # We fix the ligand. Its bond orders are correct now and it contains
    # hydrogens
    ligand.fix_bond_order(
        smiles="c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)C[C@H](CC[C@@H](C(=O)O)N)N)O)O)N",
    )
    ligand.add_hydrogens()
    assert ligand.n_atoms == 50

    # We extract the binding site using the centroid of the ligand
    binding_site = protein.extract_binding_site(ligand=ligand)

    # Extract chemical features and visualize them
    ligand_feats = oph.extract_chem_feats(ligand)
    receptor_feats = oph.extract_chem_feats(binding_site)

    viewer = oph.Viewer(binding_site=binding_site, ligands=ligand)
    viewer.add_chem_feats([ligand_feats, receptor_feats])
    viewer.show()
