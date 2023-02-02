import openpharmacophore as oph


def test_dynamic_ligand_receptor_pharmacophore(traj_er_alpha):
    # We obtain pharmacophores from a md trajectory of er-alpha
    # that consists of three frames
    protein = oph.load(traj_er_alpha)
    # We know that the file contains a single ligand
    lig_ids = protein.ligand_ids
    assert len(lig_ids) == 1
    # The receptor already contains hydrogens
    assert protein.has_hydrogens

    ligand = protein.get_ligand(lig_ids[0], remove_hyd=False)
    # Ligand already contains hydrogens
    assert ligand.has_hydrogens
    assert ligand.n_atoms == 44
    ligand.fix_bond_order(
        smiles="C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O"
    )
    assert ligand.has_aromatic_bonds()

    # We need to extract the binding site from the protein, so we can get
    # pharmacophoric features
    bsite = oph.ComplexBindingSite(protein, ligand)

    # We extract the pharmacophore
    pharmacophore = oph.LigandReceptorPharmacophore(bsite, ligand)
    pharmacophore.extract(frames=range(3))
    assert len(pharmacophore) == 3
    assert len(pharmacophore[0]) > 0
    assert len(pharmacophore[1]) > 0
    assert len(pharmacophore[2]) > 0

    # Finally we visualize the pharmacophore.
    viewer = oph.Viewer()
    viewer.add_components([protein, ligand])
    viewer.show()
