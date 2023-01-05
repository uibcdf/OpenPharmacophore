from openpharmacophore import PLComplex


def test_pl_complex_preparation(pdb_1m7w):
    # We want to prepare the protein ligand complex of the pdb 1M7W
    # to be ready for pharmacophore extraction. This complex contains
    # lauric acid (id DAO, C12H24O2) as its ligand.
    pl = PLComplex(pdb_1m7w)
    n_atoms_start = pl.topology.n_atoms
    assert pl.ligand_ids == ["DAO:B"]
    assert not pl.has_hydrogens()

    pl.set_ligand(pl.ligand_ids[0])
    # First we need to extract the ligand, so we can fix its bond orders
    # and add hydrogens to it.
    pl.ligand_to_mol()
    assert pl.ligand.GetNumAtoms() == 14
    # The ligand does not have hydrogens yet
    assert len([a for a in pl.ligand.GetAtoms() if a.GetSymbol() == "H"]) == 0
    # All bonds are single
    assert all([b.GetBondTypeAsDouble() == 1.0 for b in pl.ligand.GetBonds()])

    # After fixing ligand, it must have a double bond and hydrogens
    pl.fix_ligand(smiles="CCCCCCCCCCCC(=O)O")
    assert pl.ligand.GetNumAtoms() == 38
    assert len([a for a in pl.ligand.GetAtoms() if a.GetSymbol() == "H"]) == 24
    assert len([b.GetBondTypeAsDouble() for b in pl.ligand.GetBonds()
                if b.GetBondTypeAsDouble() == 2.0]) == 1

    # We remove the unfixed ligand from the complex
    pl.remove_ligand()
    assert pl.topology.n_atoms == n_atoms_start - 14
    n_atoms = pl.topology.n_atoms

    # We add hydrogens to the receptor
    pl.add_hydrogens()
    assert pl.has_hydrogens()
    assert pl.topology.n_atoms > n_atoms
    n_atoms = pl.topology.n_atoms

    # Finally we concatenate the fixed ligand with the receptor
    pl.add_fixed_ligand()
    assert pl.topology.n_atoms == n_atoms + 38

    pl.get_non_hyd_indices()
    assert len(pl._non_hyd_indices) == n_atoms_start

    # We create a view of our protein-ligand complex
    view = pl.show()
    assert ["nglview.adaptor.MDTrajTrajectory"] == view._ngl_component_names


def test_prepare_add_hydrogens(small_pdb_with_ligand):
    # We prepare the protein ligand complex using the method prepare of the PLComplex
    # class. In this case we want to add hydrogens to the complex.
    pl_complex = PLComplex(small_pdb_with_ligand)
    pl_complex.prepare(
        "EST:B", smiles="C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O",
        add_hydrogens=True)
    assert pl_complex.has_hydrogens()
    assert pl_complex.ligand.GetNumAtoms() == 44
    # Ligand must have aromatic bonds
    assert any(b.GetIsAromatic() for b in pl_complex.ligand.GetBonds())


def test_prepare_without_adding_hydrogens(small_pdb_with_ligand):
    # We prepare the protein ligand complex using the method prepare of the PLComplex
    # class. In this case we do not want to add hydrogens to the complex.
    pl_complex = PLComplex(small_pdb_with_ligand)
    pl_complex.prepare(
        "EST:B", smiles="C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O",
        add_hydrogens=False)
    assert not pl_complex.has_hydrogens()
    assert pl_complex.ligand.GetNumAtoms() == 20
    # Ligand must have aromatic bonds
    assert any(b.GetIsAromatic() for b in pl_complex.ligand.GetBonds())
