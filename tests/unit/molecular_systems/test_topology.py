from openpharmacophore.molecular_systems import Topology
import pytest


@pytest.fixture
def topology_2_chains():
    top = Topology()
    top.set_num_chains(2)
    top.add_atoms_to_chain(
        {
            "ALA": [("CA", "C"), ("O", "O"), ("NA", "N")],
            "MET": [("CA", "C"), ("O", "O"), ("CB", "C")],
        },
        0
    )
    top.add_atoms_to_chain(
        {
            "PRO": [("CA", "C"), ("O", "O"), ("NA", "N")],
            "LEU": [("CA", "C"), ("O", "O"), ("CB", "C")],
        },
        1
    )
    assert top.n_residues == 4
    assert top.n_atoms == 12
    return top


@pytest.fixture
def topology_with_hydrogen():
    top = Topology()
    top.set_num_chains(2)
    top.add_atoms_to_chain(
        {
            "ALA": [("CA", "C"), ("O", "O"), ("NA", "N"), ("H", "H")],
            "MET": [("CA", "C"), ("O", "O"), ("CB", "C"), ("H", "H")],
        },
        0
    )
    top.add_atoms_to_chain(
        {
            "PRO": [("CA", "C"), ("O", "O"), ("NA", "N"), ("H", "H")],
            "LEU": [("CA", "C"), ("O", "O"), ("CB", "C"), ("H", "H")],
        },
        1
    )
    assert top.n_residues == 4
    assert top.n_atoms == 16
    return top


def test_create_topology_from_scratch():
    top = Topology()
    top.set_num_chains(2)
    assert top.n_chains == 2

    res_ind_1 = top.add_residue("ALA", chain=0)
    res_ind_2 = top.add_residue("MET", chain=1)
    assert res_ind_1 == 0
    assert res_ind_2 == 1
    assert top.n_residues == 2

    atom_ind_1 = top.add_atom("CA", "C", residue=res_ind_1)
    atom_ind_2 = top.add_atom("H", "H", residue=res_ind_2)
    assert atom_ind_1 == 0
    assert atom_ind_2 == 1
    assert top.n_atoms == 2


def test_has_hydrogens(topology_2_chains,
                       topology_with_hydrogen):
    assert not topology_2_chains.has_hydrogens()
    assert topology_with_hydrogen.has_hydrogens()
