import openpharmacophore as oph


def test_pl_complex_preparation(pdb_1m7w):
    # We want to prepare the protein ligand complex of the pdb 1M7W
    # to be ready for pharmacophore extraction. This complex contains
    # lauric acid (id DAO, C12H24O2) as its ligand.
    protein = oph.load(pdb_1m7w)
    n_atoms_start = protein.n_atoms
    assert protein.has_ligands

    lig_ids = protein.ligand_ids
    assert lig_ids == ["DAO:B"]
    assert not protein.has_hydrogens

    # First we need to extract the ligand, so we can fix its bond orders
    # and add hydrogens to it.
    ligand = protein.get_ligand(lig_ids[0], remove=True)
    assert ligand.n_atoms == 14
    # The protein should have 14 less atoms because the ligand was removed
    assert protein.n_atoms == n_atoms_start - 14
    # The ligand does not have hydrogens yet and all its bonds are single
    assert not ligand.has_hydrogens()
    assert not ligand.has_double_bonds()

    # After fixing ligand, it must have a double bond and hydrogens
    ligand.fix_bond_order(smiles="CCCCCCCCCCCC(=O)O")
    assert ligand.has_double_bonds()

    ligand.add_hydrogens()
    assert ligand.n_atoms == 38
    assert ligand.has_hydrogens()

    # We add hydrogens to the receptor
    protein.add_hydrogens()
    assert protein.has_hydrogens()
